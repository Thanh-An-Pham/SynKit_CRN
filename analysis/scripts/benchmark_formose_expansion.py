"""Run isolated, sequential expansion benchmarks with process-group limits."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import signal
import statistics
import subprocess
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--synkit-python', required=True)
    p.add_argument('--mod-python', required=True)
    p.add_argument('--molfoundry-python', required=True)
    p.add_argument('--synkit-root',required=True)
    p.add_argument('--expander',required=True)
    p.add_argument('--depths',type=int,nargs='+',default=list(range(1,8)))
    p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--timeout',type=float,default=180)
    p.add_argument('--memory-gib',type=int,default=8)
    p.add_argument('--output-dir',type=Path,default=Path('analysis/formose-expansion'))
    args=p.parse_args()
    out=args.output_dir.resolve();(out/'runs').mkdir(parents=True,exist_ok=True)
    worker=Path(__file__).with_name('formose_expansion_worker.py').resolve()
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
             PYTHONPATH=args.synkit_root,PYTHONHASHSEED='0')
    cpu=min(os.sched_getaffinity(0))
    def limits():
        os.sched_setaffinity(0,{cpu})
        cap=args.memory_gib*1024**3
        resource.setrlimit(resource.RLIMIT_AS,(cap,cap))
    report={'schema':'synkit.formose-expansion-comparison/1',
            'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'platform':platform.platform(),'cpu':next(x.split(':',1)[1].strip() for x in Path('/proc/cpuinfo').read_text().splitlines() if x.startswith('model name')),
            'cpu_affinity':[cpu],'threads':1,'timeout_seconds_per_worker':args.timeout,
            'memory_gib_per_worker':args.memory_gib,'repeats_requested':args.repeats,
            'warmups_per_worker':1,'arguments':{k:str(v) for k,v in vars(args).items()},
            'timing_scope':'Cumulative generation from seeds, including rule/seed parsing and network construction; excludes imports, engine warm-up, GML translation, canonicalization audit and structural analyses.',
            'policy':'Both seeds initially active; at least one frontier reactant per iteration; no repeated species in a reactant mixture; no chemical size cutoff; SynKit task/mixture caps raised to 2**60.',
            'source_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [worker,Path(__file__),out/'inputs.json']},'runs':[]}
    for depth in args.depths:
        for tool in ['synkit','mod','molfoundry']:
            dest=out/'runs'/f'{tool}-{depth}.json'; log=dest.with_suffix('.log')
            cmd=[getattr(args,tool+'_python'),str(worker),'--tool',tool,'--depth',str(depth),
                 '--repeats',str(args.repeats),'--inputs',str(out/'inputs.json'),'--output',str(dest),
                 '--synkit-root',args.synkit_root,'--expander',args.expander]
            print(f'Start {tool} depth {depth}',flush=True)
            started=time.perf_counter();status='complete'
            with log.open('w') as stream:
                proc=subprocess.Popen(cmd,stdout=stream,stderr=subprocess.STDOUT,env=env,
                                      start_new_session=True,preexec_fn=limits)
                try:
                    code=proc.wait(timeout=args.timeout)
                    if code: status='error'
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid,signal.SIGKILL);proc.wait();status='timeout';code=proc.returncode
            row={'tool':tool,'depth':depth,'status':status,'returncode':code,'worker_wall_seconds':time.perf_counter()-started,
                 'result':str(dest.relative_to(out)) if dest.exists() else None,'log':str(log.relative_to(out))}
            if dest.exists():
                data=json.loads(dest.read_text()); samples=data['samples']
                row.update(n_species=data['n_species'],n_reactions=data['n_reactions'],samples=samples,
                           snapshot_sha256=data['snapshot']['sha256'],version=data['version'])
                if samples:
                    q=statistics.quantiles(samples,n=4,method='inclusive') if len(samples)>1 else [samples[0]]*3
                    row.update(median_seconds=statistics.median(samples),q1_seconds=q[0],q3_seconds=q[2])
            report['runs'].append(row)
            peers=[r for r in report['runs'] if r['depth']==depth and r['status']=='complete']
            report.setdefault('equivalence',{})[str(depth)]={
                'tools_complete':[r['tool'] for r in peers],
                'all_three_complete':len(peers)==3,
                'identical_snapshots':len(peers)>1 and len({r['snapshot_sha256'] for r in peers})==1}
            (out/'comparison.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
            print(json.dumps(row),flush=True)
    print('Wrote '+str(out/'comparison.json'),flush=True)


if __name__=='__main__':
    main()
