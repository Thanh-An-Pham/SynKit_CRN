"""Sequential SynKit variant controls, preserving the original tool comparison."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import statistics
import subprocess
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--depths',type=int,nargs='+',default=list(range(1,8)))
    p.add_argument('--timeout',type=float,default=180)
    p.add_argument('--wait-pid',type=int)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[2]/'analysis/formose-expansion'
    original=json.loads((root/'comparison.json').read_text());config=original['arguments']
    worker=Path(__file__).with_name('formose_synkit_variants.py')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONHASHSEED='0',PYTHONPATH=config['synkit_root'])
    cpu=original['cpu_affinity'][0]
    def limits():
        os.sched_setaffinity(0,{cpu});cap=8*1024**3
        resource.setrlimit(resource.RLIMIT_AS,(cap,cap))
    if a.wait_pid:
        print('Waiting for the original benchmark to finish before timing controls.',flush=True)
        while Path(f'/proc/{a.wait_pid}').exists(): time.sleep(1)
    report={'schema':'synkit.formose-variants/1','cpu_affinity':[cpu],'timeout_seconds_per_worker':a.timeout,
            'source_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),worker]},'runs':[]}
    for depth in a.depths:
        # Include native Lewis controls again, to measure the native-request split.
        for mode in ['python-explicit','native-explicit','native-lewis']:
            if mode=='native-lewis' and depth!=6: continue
            dest=root/'runs'/f'{mode}-{depth}.json';log=dest.with_suffix('.log')
            cmd=[config['synkit_python'],str(worker),'--mode',mode,'--depth',str(depth),'--inputs',str(root/'inputs.json'),
                 '--output',str(dest),'--synkit-root',config['synkit_root'],'--expander',config['expander']]
            print(f'Start {mode} depth {depth}',flush=True);started=time.perf_counter()
            status='complete'
            with log.open('w') as stream:
                proc=subprocess.Popen(cmd,stdout=stream,stderr=subprocess.STDOUT,env=env,start_new_session=True,preexec_fn=limits)
                try:
                    code=proc.wait(timeout=a.timeout)
                    if code:status='error'
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid,signal.SIGKILL);proc.wait();status='timeout';code=proc.returncode
            row={'mode':mode,'depth':depth,'status':status,'returncode':code,'worker_wall_seconds':time.perf_counter()-started,
                 'result':str(dest.relative_to(root)) if dest.exists() else None,'log':str(log.relative_to(root))}
            if dest.exists():
                data=json.loads(dest.read_text());samples=data['samples']
                row.update(samples=samples,n_species=data['n_species'],n_reactions=data['n_reactions'],snapshot_sha256=data['snapshot']['sha256'])
                peer=next((r for r in original['runs'] if r['tool']=='synkit' and r['depth']==depth and r['status']=='complete'),None)
                row['matches_original_native_lewis']=peer is not None and row['snapshot_sha256']==peer['snapshot_sha256']
                if peer:assert row['matches_original_native_lewis'],(mode,depth,'network mismatch')
                if samples:
                    q=statistics.quantiles(samples,n=4,method='inclusive') if len(samples)>1 else [samples[0]]*3
                    row.update(median_seconds=statistics.median(samples),q1_seconds=q[0],q3_seconds=q[2])
                if mode.startswith('native') and samples:
                    row['median_native_request_seconds']=statistics.median(s['native_request_including_transport_seconds'] for s in data['component_samples'])
                    row['median_python_construction_seconds']=statistics.median(s['python_construction_seconds'] for s in data['component_samples'])
            report['runs'].append(row)
            (root/'variants.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
            print(json.dumps(row),flush=True)


if __name__=='__main__':main()
