"""Hydrogen/implementation controls for the matched formose expansion."""
import argparse
from collections import Counter
import contextlib
import json
import logging
from pathlib import Path
import sys
import time

from formose_expansion_worker import normalize


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['python-explicit','native-explicit','native-lewis'],required=True)
    p.add_argument('--depth',type=int,required=True)
    p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--inputs',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--synkit-root',required=True)
    p.add_argument('--expander',required=True)
    a=p.parse_args();sys.path.insert(0,a.synkit_root)
    from synkit.CRN import CRNExpand
    from synkit.CRN.Construct import native_expand as ne
    import synkit
    from rdkit import rdBase
    logging.disable(logging.INFO)
    inputs=json.loads(Path(a.inputs).read_text())
    native=a.mode.startswith('native')
    hydrogen='lewis' if a.mode.endswith('lewis') else 'explicit'
    result={'mode':a.mode,'depth':a.depth,'synkit_version':synkit.__version__,
            'rdkit_python_version':rdBase.rdkitVersion,'samples':[],'component_samples':[]}
    original=ne.run_native_expansion
    request_seconds=[]
    def measured(*args,**kwargs):
        started=time.perf_counter();value=original(*args,**kwargs)
        request_seconds.append(time.perf_counter()-started)
        return value
    ne.run_native_expansion=measured
    context=ne.native_expander(a.expander) if native else contextlib.nullcontext()
    with context:
        for repeat in range(a.repeats+1):
            request_seconds.clear()
            started=time.perf_counter()
            builder=CRNExpand(rules=inputs['rules'],repeats=a.depth,keep_aam=False,
                              rule_backend='native' if native else 'rdkit',hydrogen_model=hydrogen,
                              threads=1,max_mixtures_per_rule_step=2**60,max_tasks_per_step=2**60)
            graph=builder.build(inputs['seeds'])
            elapsed=time.perf_counter()-started
            print(json.dumps({'phase':'generation_complete','mode':a.mode,'depth':a.depth,'repeat':repeat,
                              'seconds':elapsed,'native_request_seconds':sum(request_seconds) if native else None}),flush=True)
            species=[d['smiles'] for _,d in graph.nodes(data=True) if d.get('kind')=='species']
            reactions=[]
            for n,d in graph.nodes(data=True):
                if d.get('kind')!='rule': continue
                left=[];right=[]
                for u,_,e in graph.in_edges(n,data=True): left.extend([graph.nodes[u]['smiles']]*int(e.get('stoich',1)))
                for _,v,e in graph.out_edges(n,data=True): right.extend([graph.nodes[v]['smiles']]*int(e.get('stoich',1)))
                reactions.append((left,right))
            snapshot=normalize(species,reactions)
            if repeat:
                assert snapshot==result['snapshot']
                result['samples'].append(elapsed)
                result['component_samples'].append({'total_seconds':elapsed,
                    'native_request_including_transport_seconds':sum(request_seconds) if native else None,
                    'python_construction_seconds':elapsed-sum(request_seconds) if native else None})
            result.update(snapshot=snapshot,n_species=len(snapshot['species']),n_reactions=len(snapshot['reactions']),
                          native_report=builder.native_report,status='running')
            if native:
                assert all(not s['capped_rules'] and not s['task_budget_exhausted'] for s in builder.native_report['steps'])
            Path(a.output).write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    result['status']='complete'
    Path(a.output).write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')


if __name__=='__main__': main()
