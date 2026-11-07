"""Compare full-API native SynReactor DPO, MØD and MolFoundry expansion."""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import resource
import select
import signal
import statistics
import subprocess
import sys
import time

from formose_expansion_worker import normalize, rule_gml

ROOT = Path(__file__).resolve().parents[2] / 'analysis/formose-expansion'


def save(path, value):
    path = Path(path)
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')
    tmp.replace(path)


def emit(value):
    print(json.dumps(value), flush=True)


def serve(a):
    """Retain imported engines, but construct every timed network from scratch."""
    logging.disable(logging.INFO)
    inputs = json.loads((ROOT / 'inputs.json').read_text())
    result = {'backend': a.backend, 'depth': a.depth, 'samples': [],
              'warmup_samples': [], 'audits': [], 'status': 'running'}
    context = contextlib.nullcontext()
    if a.backend == 'synreactor-native':
        sys.path.insert(0, a.synkit_root)
        from synkit.CRN import CRNExpand
        import synkit
        result.update(version=synkit.__version__, radical_policy='ignore',
                      controller='cpp', kernel='cpp',
                      hydrogen_representation='implicit counts')
        from synkit.CRN.Construct.native_dpo import native_dpo_expander
        context = native_dpo_expander(a.dpo_expander)

        def expand():
            builder = CRNExpand(
                rules=inputs['rules'], repeats=a.depth, rule_backend='synreactor-native',
                native_dpo_expander=a.dpo_expander,
                dpo_radical_policy='ignore', keep_aam=False, explicit_h=False,
                implicit_temp=True, hydrogen_model='lewis', threads=1,
                max_mixtures_per_rule_step=2**60, max_tasks_per_step=2**60)
            return builder, builder.build(inputs['seeds'])

        def export(value):
            builder, graph = value
            species = [d['smiles'] for _, d in graph.nodes(data=True) if d.get('kind') == 'species']
            reactions = []
            for n, d in graph.nodes(data=True):
                if d.get('kind') != 'rule':
                    continue
                left, right = [], []
                for u, _, e in graph.in_edges(n, data=True):
                    left.extend([graph.nodes[u]['smiles']] * int(e.get('stoich', 1)))
                for _, v, e in graph.out_edges(n, data=True):
                    right.extend([graph.nodes[v]['smiles']] * int(e.get('stoich', 1)))
                reactions.append((left, right))
            result['dpo_report'] = builder.dpo_report
            return species, reactions
    else:
        if a.backend == 'mod':
            import mod as api
            result['version'] = str(api.version())
        else:
            import molfoundry as api
            result['version'] = importlib.metadata.version('molfoundry')
        result['hydrogen_representation'] = 'explicit vertices'
        # Fixed format translation is outside the API timing; rule parsing is inside.
        gml = [rule_gml(s, f'r{i}') for i, s in enumerate(inputs['rules'])]

        def expand():
            if a.backend == 'mod':
                seeds = [api.Graph.fromSMILES(s, add=False) for s in inputs['seeds']]
                rules = [api.Rule.fromGMLString(s, add=False) for s in gml]
            else:
                seeds = [api.smiles(s, add=False) for s in inputs['seeds']]
                rules = [api.ruleGMLString(s, add=False) for s in gml]
            strategy = api.addSubset(seeds) >> api.repeat[a.depth](api.leftPredicate[
                lambda d: len(set(g.id for g in d.left)) == len(d.left)
            ](rules))
            dg = api.DG(graphDatabase=seeds)
            with dg.build() as builder:
                builder.execute(strategy, verbosity=0)
            return dg

        def export(dg):
            return ([v.graph.smiles for v in dg.vertices],
                    [([v.graph.smiles for v in e.sources], [v.graph.smiles for v in e.targets])
                     for e in dg.edges])

    expected_path = ROOT / 'runs' / f'synkit-{a.depth}.json'
    expected = json.loads(expected_path.read_text())['snapshot'] if expected_path.exists() else None
    with context:
        emit({'phase': 'ready', 'version': result['version']})
        for line in sys.stdin:
            command = json.loads(line)
            if command['phase'] == 'stop':
                break
            # Tool chatter is kept in the diagnostic log, away from the protocol.
            with contextlib.redirect_stdout(sys.stderr):
                started = time.perf_counter()
                value = expand()
                elapsed = time.perf_counter() - started
            emit({'phase': 'generation_complete', 'seconds': elapsed})
            with contextlib.redirect_stdout(sys.stderr):
                snapshot = normalize(*export(value))
            del value
            if expected is not None:
                assert snapshot == expected, 'Completed network differs from matched reference'
            if 'snapshot' in result:
                assert snapshot == result['snapshot'], 'Network changed between repetitions'
            result.update(snapshot=snapshot, n_species=len(snapshot['species']),
                          n_reactions=len(snapshot['reactions']))
            result['warmup_samples' if command['phase'] == 'warmup' else 'samples'].append(elapsed)
            result['audits'].append({'phase': command['phase'], 'index': command['index'],
                                     'snapshot_sha256': snapshot['sha256']})
            save(a.output, result)
            emit({'phase': 'done', 'seconds': elapsed, 'snapshot_sha256': snapshot['sha256']})
    result['status'] = 'complete'
    save(a.output, result)


class Worker:
    def __init__(self, command, log, env, limits):
        self.log = log.open('w')
        self.proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=self.log, env=env, start_new_session=True,
                                     preexec_fn=limits, bufsize=0)
        self.buffer = b''

    def receive(self, timeout):
        deadline = time.monotonic() + timeout
        while b'\n' not in self.buffer:
            remaining = deadline - time.monotonic()
            if remaining <= 0 or not select.select([self.proc.stdout], [], [], remaining)[0]:
                raise TimeoutError
            chunk = os.read(self.proc.stdout.fileno(), 65536)
            if not chunk:
                raise RuntimeError(f'Worker exited: {self.proc.poll()}')
            self.buffer += chunk
        line, self.buffer = self.buffer.split(b'\n', 1)
        message = json.loads(line)
        self.log.write(json.dumps(message) + '\n')
        self.log.flush()
        return message

    def send(self, command):
        self.proc.stdin.write((json.dumps(command) + '\n').encode())
        self.proc.stdin.flush()

    def kill(self):
        try:
            os.killpg(self.proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        self.proc.wait()
        self.log.close()

    def finish(self):
        self.send({'phase': 'stop'})
        code = self.proc.wait(timeout=30)
        self.log.close()
        assert code == 0, f'Worker exited: {code}'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--synkit-root', required=True)
    p.add_argument('--library', required=True)
    p.add_argument('--dpo-expander', required=True)
    p.add_argument('--external', action='store_true')
    p.add_argument('--mod-python', default='/tmp/formose-mod-env/bin/python')
    p.add_argument('--molfoundry-python', default='/tmp/formose-comparison-env/bin/python')
    p.add_argument('--depths', type=int, nargs='+', default=list(range(1, 7)))
    p.add_argument('--repeats', type=int, default=20)
    p.add_argument('--warmups', type=int, default=3)
    p.add_argument('--timeout', type=float, default=180)
    p.add_argument('--worker', action='store_true')
    p.add_argument('--backend', choices=['synreactor-native', 'mod', 'molfoundry'])
    p.add_argument('--depth', type=int)
    p.add_argument('--output')
    a = p.parse_args()
    if a.worker:
        return serve(a)
    assert a.repeats >= 2 and a.warmups >= 1
    out = ROOT / 'dpo'
    out.mkdir(exist_ok=True)
    cpu = min(os.sched_getaffinity(0))

    def limits():
        os.sched_setaffinity(0, {cpu})
        cap = 8 * 1024**3
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))

    def source_hashes():
        checkout = Path(a.synkit_root)
        return {str(f.relative_to(checkout)): hashlib.sha256(f.read_bytes()).hexdigest()
                for directory in [checkout/'synkit', checkout/'Experiment/CRN']
                for f in sorted(directory.rglob('*'))
                if f.is_file() and f.suffix in {'.py', '.cpp', '.h', '.hpp'}}

    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               MKL_NUM_THREADS='1', PYTHONHASHSEED='0', PYTHONPATH=a.synkit_root)
    sys.path.insert(0, a.synkit_root)
    from synkit.CRN.Construct.native_dpo import dpo_expander_manifest
    script = Path(__file__).read_bytes()
    digest = hashlib.sha256(script).hexdigest()
    script_source = out/'provenance/drivers'/f'{digest}.py'
    script_source.parent.mkdir(parents=True, exist_ok=True)
    script_source.write_bytes(script)
    report = {'schema': 'synkit.formose-dpo-comparison/4', 'threads': 1,
              'cpu_affinity': [cpu], 'timeout_seconds_per_expansion': a.timeout,
              'memory_gib_per_process': 8, 'repeats_requested': a.repeats,
              'warmups': a.warmups, 'source_sha256': source_hashes(),
              'native_library': a.library,
              'native_manifest': json.loads(Path(a.library).with_suffix('.json').read_text()),
              'native_dpo_expander': a.dpo_expander,
              'native_dpo_manifest': dpo_expander_manifest(a.dpo_expander),
              'script_sha256': digest, 'script_source': str(script_source.relative_to(ROOT)),
              'timing_scope': 'Fresh rule/seed preparation, matching and DPO rewriting, molecular identity and network construction through each tool API. Imports, library startup, fixed format translation and independent identity audits excluded. Native controller retained across calls.',
              'order_policy': 'Only one active expansion. Cyclic rotation across tools in every warm-up and measured round; starting tool also rotates across depths.',
              'radical_policy': 'ignore: generic formose patterns omit radical guards; kernel retains zero radical deltas',
              'schedule': [], 'runs': []}
    save(out/'comparison.json', report)
    for depth in a.depths:
        names = ['synreactor-native'] + (['mod', 'molfoundry'] if a.external else [])
        workers, rows = {}, {}
        try:
            for backend in names:
                dest = out/f'{backend}-{depth}.json'
                dest.unlink(missing_ok=True)
                interpreter = getattr(a, backend+'_python') if backend in {'mod', 'molfoundry'} else sys.executable
                command = [interpreter, str(Path(__file__).resolve()), '--worker',
                           '--synkit-root', a.synkit_root, '--library', a.library,
                           '--dpo-expander', a.dpo_expander, '--backend', backend,
                           '--depth', str(depth), '--output', str(dest)]
                row = {'backend': backend, 'depth': depth, 'status': 'running',
                       'result': dest.name, 'log': dest.with_suffix('.log').name}
                rows[backend] = row
                report['runs'].append(row)
                w = Worker(command, dest.with_suffix('.log'), env, limits)
                workers[backend] = w
                message = w.receive(180)
                assert message['phase'] == 'ready'
                row['version'] = message['version']
            for round_index in range(a.warmups+a.repeats):
                phase = 'warmup' if round_index < a.warmups else 'measured'
                index = round_index if phase == 'warmup' else round_index-a.warmups
                shift = (round_index+depth-1) % len(names)
                for backend in names[shift:]+names[:shift]:
                    if backend not in workers:
                        continue
                    w, row = workers[backend], rows[backend]
                    event = {'depth': depth, 'backend': backend, 'phase': phase, 'index': index}
                    report['schedule'].append(event)
                    started = time.monotonic()
                    try:
                        w.send({'phase': phase, 'index': index})
                        message = w.receive(a.timeout)
                        assert message['phase'] == 'generation_complete'
                        event['seconds'] = message['seconds']
                        # Separate budget: audits never consume an expansion's limit.
                        audited = w.receive(180)
                        assert audited['phase'] == 'done'
                        event['snapshot_sha256'] = audited['snapshot_sha256']
                    except (TimeoutError, RuntimeError, AssertionError) as error:
                        row.update(status='timeout' if isinstance(error, TimeoutError) else 'error',
                                   failed_phase=phase, failed_index=index,
                                   failure_stage='audit' if 'seconds' in event else 'expansion',
                                   failure=str(error) or 'Per-expansion time limit reached')
                        w.kill()
                        del workers[backend]
                        print(json.dumps(row), flush=True)
                    event['wall_seconds_with_audit'] = time.monotonic()-started
                    if (out/row['result']).exists():
                        data = json.loads((out/row['result']).read_text())
                        samples = data['samples']
                        row.update(samples=samples, warmup_samples=data['warmup_samples'],
                                   n_species=data['n_species'], n_reactions=data['n_reactions'],
                                   snapshot_sha256=data['snapshot']['sha256'])
                        if samples:
                            q1, _, q3 = statistics.quantiles(samples, n=4, method='inclusive') if len(samples)>1 else [samples[0]]*3
                            row.update(median_seconds=statistics.median(samples), q1_seconds=q1, q3_seconds=q3)
                    save(out/'comparison.json', report)
                print(f'Depth {depth}: {phase} {index+1} complete', flush=True)
            for backend, w in list(workers.items()):
                w.finish()
                del workers[backend]
                rows[backend]['status'] = 'complete'
                print(json.dumps(rows[backend]), flush=True)
        finally:
            for w in workers.values():
                w.kill()
        save(out/'comparison.json', report)
    assert report['source_sha256'] == source_hashes(), 'Sources changed during comparison'
    report['source_unchanged'] = True
    save(out/'comparison.json', report)


if __name__ == '__main__':
    main()
