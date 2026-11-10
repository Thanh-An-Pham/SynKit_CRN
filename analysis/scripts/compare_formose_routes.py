"""Independently check the ten prescribed formose routes in MolFoundry."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
from collections import Counter

import molfoundry as mf
import os
import sys
from rdkit import Chem
sys.path.insert(0, str(Path(os.environ.get('SYNKIT_ROOT', Path(__file__).resolve().parents[4] / 'SynEco/SynKit'))))
from Experiment.CRN.studies.formose_dpo import rule_gml

def canonical(s):
    from rdkit import Chem
    m = Chem.MolFromSmiles(s)
    if m is None:
        raise ValueError(f'Invalid exported molecule: {s}')
    return Chem.MolToSmiles(Chem.RemoveHs(m), canonical=True, isomericSmiles=False)


def main():
    paper = Path(__file__).resolve().parents[2]
    root = paper/'analysis/formose-expansion'
    inputs = json.loads((root/'dpo/inputs.json').read_text())
    gs = [mf.smiles(s, add=False) for s in inputs['seeds']]
    rules = [mf.ruleGMLString(rule_gml(s,f'r{i}'),add=False) for i,s in enumerate(inputs['rules'])]
    dg = mf.DG(graphDatabase=gs)
    with dg.build() as builder:
        builder.execute(mf.addSubset(gs) >> mf.repeat[4](mf.leftPredicate[
            lambda d: len(set(g.id for g in d.left)) == len(d.left)
        ](rules)),verbosity=0)
    vertices = {canonical(v.graph.smiles):v for v in dg.vertices}
    c2,enol,t4,t6 = map(canonical,['O=CCO','OC=CO','O=CC(O)C(O)CO','O=CC(O)C(O)C(O)C(O)CO'])
    definitions = {'r_1':([c2],[enol]), 'r_3':([enol,c2],[t4]), 'r_7':([enol,t4],[t6])}
    edges = {}
    for name,(left,right) in definitions.items():
        found = [e for e in dg.edges if Counter(canonical(v.graph.smiles) for v in e.sources)==Counter(left)
                 and Counter(canonical(v.graph.smiles) for v in e.targets)==Counter(right)]
        assert len(found)==1,(name,len(found))
        edges[name]=found[0]
    rows=[]
    original=paper/'analysis/crn-supplement-audit.json'
    for row in json.loads(original.read_text())['formose_routes']:
        outcome=row['outcome']; pool=outcome['initial_marking']; flux=outcome['requested_flux']
        solution=mf.FlowSolution(dg,[(edges[n],c) for n,c in flux.items()],
                                 [(vertices[canonical(s)],c) for s,c in pool.items()],[],0)
        answer=solution.is_realizable()
        cert=next(solution.enumerate_certificates(limit=1),None)
        assert (cert is not None)==answer
        replay=None;sequence=None
        if cert is not None:
            marking=Counter({canonical(s):c for s,c in pool.items()})
            reverse={e.id:n for n,e in edges.items()};sequence=[]
            for e in cert.edges:
                left=Counter(canonical(v.graph.smiles) for v in e.sources)
                right=Counter(canonical(v.graph.smiles) for v in e.targets)
                assert all(marking[s]>=c for s,c in left.items())
                marking.subtract(left);marking.update(right);sequence.append(reverse[e.id])
            assert Counter(sequence)==Counter(flux)
            assert marking[t4 if row['target']=='T4' else t6]>=1
            replay=True
        assert answer==outcome['realizable']
        rows.append({'pool':row['pool'],'target':row['target'],'initial_marking':pool,
                     'requested_flux':flux,'synkit_realizable':outcome['realizable'],
                     'molfoundry_realizable':answer,'certificate':sequence,'independent_replay_ok':replay})
    result={'schema':'synkit.formose-route-comparison/1','molfoundry_version':importlib.metadata.version('molfoundry'),
            'all_ten_agree':len(rows)==10,'positive_certificates_replayed':sum(r['independent_replay_ok'] is True for r in rows),
            'scope':'Fixed reaction multisets from finite initial markings; does not test unrestricted target reachability. FlowSolution is constructed directly, including infeasible candidates, without solving an aggregate hyperflow model.',
            'rows':rows}
    (root/'route-comparison.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
