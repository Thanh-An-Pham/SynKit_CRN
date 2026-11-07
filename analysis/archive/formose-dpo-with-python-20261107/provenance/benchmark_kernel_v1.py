"""Benchmark SynReactor's Python and C++ DPO kernels with a Python controller."""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
import resource
import signal
import statistics
import subprocess
import sys
import time

from formose_expansion_worker import normalize


def worker(a):
    sys.path.insert(0,a.synkit_root)
    from synkit.CRN import CRNExpand
    logging.disable(logging.INFO)
    root=Path(__file__).resolve().parents[2]/'analysis/formose-expansion'
    inputs=json.loads((root/'inputs.json').read_text())
    result={'backend':a.backend,'depth':a.depth,'radical_policy':'ignore',
            'controller':'python','kernel':'cpp' if a.backend=='synreactor-cpp' else 'python',
            'engine':'SynReactor DPO; full matches; implicit H', 'samples':[], 'status':'running'}
    for repeat in range(a.repeats+1):
        started=time.perf_counter()
        builder=CRNExpand(rules=inputs['rules'],repeats=a.depth,rule_backend=a.backend,
                          reactor_library=a.library if a.backend=='synreactor-cpp' else None,
                          dpo_radical_policy='ignore',keep_aam=False,explicit_h=False,
                          implicit_temp=True,hydrogen_model='lewis',threads=1,
                          max_mixtures_per_rule_step=2**60,max_tasks_per_step=2**60)
        graph=builder.build(inputs['seeds'])
        elapsed=time.perf_counter()-started
        print(json.dumps({'phase':'generation_complete','backend':a.backend,'depth':a.depth,
                          'repeat':repeat,'seconds':elapsed}),flush=True)
        species=[d['smiles'] for _,d in graph.nodes(data=True) if d.get('kind')=='species']
        reactions=[]
        for n,d in graph.nodes(data=True):
            if d.get('kind')!='rule':continue
            left=[];right=[]
            for u,_,e in graph.in_edges(n,data=True):left.extend([graph.nodes[u]['smiles']]*int(e.get('stoich',1)))
            for _,v,e in graph.out_edges(n,data=True):right.extend([graph.nodes[v]['smiles']]*int(e.get('stoich',1)))
            reactions.append((left,right))
        snap=normalize(species,reactions)
        if repeat:
            assert snap==result['snapshot']
            result['samples'].append(elapsed)
        previous=root/'runs'/f'synkit-{a.depth}.json'
        if previous.exists():
            assert snap==json.loads(previous.read_text())['snapshot'],'DPO network differs from matched formose benchmark'
        result.update(snapshot=snap,n_species=len(snap['species']),n_reactions=len(snap['reactions']),dpo_report=builder.dpo_report)
        Path(a.output).write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
        print(json.dumps({'phase':'network_verified','species':result['n_species'],'reactions':result['n_reactions']}),flush=True)
    result['status']='complete'
    Path(a.output).write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--synkit-root',required=True)
    p.add_argument('--library',required=True)
    p.add_argument('--depths',type=int,nargs='+',default=list(range(1,8)))
    p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--timeout',type=float,default=180)
    p.add_argument('--worker',action='store_true')
    p.add_argument('--backend',choices=['synreactor-python','synreactor-cpp'])
    p.add_argument('--depth',type=int)
    p.add_argument('--output')
    a=p.parse_args()
    if a.worker:return worker(a)
    root=Path(__file__).resolve().parents[2]/'analysis/formose-expansion'
    out=root/'dpo';out.mkdir(exist_ok=True)
    cpu=min(os.sched_getaffinity(0))
    def limits():
        os.sched_setaffinity(0,{cpu});cap=8*1024**3
        resource.setrlimit(resource.RLIMIT_AS,(cap,cap))
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONHASHSEED='0',PYTHONPATH=a.synkit_root)
    report={'schema':'synkit.formose-dpo-comparison/1','controller':'python',
            'threads':1,'cpu_affinity':[cpu],'timeout_seconds_per_worker':a.timeout,
            'memory_gib_per_process':8,'repeats_requested':a.repeats,'warmups':1,
            'radical_policy':'ignore: generic formose patterns omit radical guards; kernel retains zero radical deltas',
            'native_library':a.library,'native_manifest':json.loads(Path(a.library).with_suffix('.json').read_text()),
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'timing_scope':'Fresh rule compilation and generator/controller construction, rule matching, DPO rewriting, output conversion and CRN construction; imports and kernel library warm-up excluded.',
            'runs':[]}
    for depth in a.depths:
        for backend in ['synreactor-python','synreactor-cpp']:
            dest=out/f'{backend}-{depth}.json';log=dest.with_suffix('.log')
            cmd=[sys.executable,str(Path(__file__).resolve()),'--worker','--synkit-root',a.synkit_root,
                 '--library',a.library,'--depth',str(depth),'--backend',backend,'--repeats',str(a.repeats),'--output',str(dest)]
            print(f'Start {backend} depth {depth}',flush=True);started=time.perf_counter();status='complete'
            with log.open('w') as stream:
                proc=subprocess.Popen(cmd,stdout=stream,stderr=subprocess.STDOUT,env=env,start_new_session=True,preexec_fn=limits)
                try:
                    code=proc.wait(timeout=a.timeout)
                    if code:status='error'
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid,signal.SIGKILL);proc.wait();code=proc.returncode;status='timeout'
            row={'backend':backend,'depth':depth,'status':status,'returncode':code,
                 'wall_seconds':time.perf_counter()-started,'result':dest.name if dest.exists() else None,'log':log.name}
            if dest.exists():
                data=json.loads(dest.read_text());samples=data['samples']
                row.update(samples=samples,n_species=data['n_species'],n_reactions=data['n_reactions'],snapshot_sha256=data['snapshot']['sha256'])
                if samples:
                    q=statistics.quantiles(samples,n=4,method='inclusive') if len(samples)>1 else [samples[0]]*3
                    row.update(median_seconds=statistics.median(samples),q1_seconds=q[0],q3_seconds=q[2])
            report['runs'].append(row)
            (out/'comparison.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
            print(json.dumps(row),flush=True)


if __name__=='__main__':main()
