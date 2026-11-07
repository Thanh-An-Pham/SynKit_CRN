"""One isolated formose expansion; timings exclude imports and identity audit."""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import sys
import time
from collections import Counter


def rule_gml(smarts, name):
    """Translate this element/bond-only mapped rule to an explicit-H DPO rule."""
    from rdkit import Chem
    sides = []
    nodes = {}
    for side in smarts.split('>>'):
        mol = Chem.MolFromSmarts(side)
        nodes.update({a.GetAtomMapNum(): a.GetSymbol() for a in mol.GetAtoms()})
        edges = {}
        for b in mol.GetBonds():
            ends = tuple(sorted((b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum())))
            edges[ends] = {1: '-', 2: '=', 3: '#'}[int(b.GetBondTypeAsDouble())]
        sides.append(edges)
    shared = {k: v for k, v in sides[0].items() if sides[1].get(k) == v}
    def edge_text(edges):
        return '\n'.join(f'edge [ source {a} target {b} label "{v}" ]' for (a,b),v in edges.items())
    return ('rule [ ruleID "' + name + '" left [\n' + edge_text({k:v for k,v in sides[0].items() if k not in shared})
            + '\n] context [\n' + '\n'.join(f'node [ id {i} label "{v}" ]' for i,v in nodes.items())
            + '\n' + edge_text(shared) + '\n] right [\n' + edge_text({k:v for k,v in sides[1].items() if k not in shared}) + '\n] ]')


def canonical(s):
    from rdkit import Chem
    m = Chem.MolFromSmiles(s)
    if m is None:
        raise ValueError(f'Invalid exported molecule: {s}')
    return Chem.MolToSmiles(Chem.RemoveHs(m), canonical=True, isomericSmiles=False)


def normalize(species, reactions):
    sp = sorted(set(map(canonical, species)))
    def side(values):
        return sorted(Counter(map(canonical, values)).items())
    rs = sorted({json.dumps([side(a), side(b)], separators=(',', ':')) for a,b in reactions})
    payload = {'species': sp, 'reactions': [json.loads(x) for x in rs]}
    payload['sha256'] = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return payload


def run(args):
    inputs = json.loads(Path(args.inputs).read_text())
    rules, seeds = inputs['rules'], inputs['seeds']
    logging.disable(logging.INFO)
    result = {'tool': args.tool, 'depth': args.depth, 'repeats': args.repeats, 'samples': []}
    if args.tool == 'synkit':
        sys.path.insert(0, args.synkit_root)
        from synkit.CRN import CRNExpand
        from synkit.CRN.Construct.native_expand import native_expander, expander_manifest
        import synkit
        result['version'] = synkit.__version__
        result['native_manifest'] = expander_manifest(args.expander)
        context = native_expander(args.expander)
    else:
        if args.tool == 'mod':
            import mod as api
            result['version'] = str(api.version())
        else:
            import molfoundry as api
            result['version'] = importlib.metadata.version('molfoundry')
        gml_rules = [rule_gml(s, f'r{i}') for i,s in enumerate(rules)]
        context = contextlib.nullcontext()

    with context:
        # One untimed warm-up: start embedded engines and the native server.
        for repeat in range(args.repeats + 1):
            started = time.perf_counter()
            if args.tool == 'synkit':
                builder = CRNExpand(rules=rules, repeats=args.depth, keep_aam=False,
                                    rule_backend='native', hydrogen_model='lewis', threads=1,
                                    max_mixtures_per_rule_step=2**60, max_tasks_per_step=2**60)
                graph = builder.build(seeds)
            else:
                if args.tool == 'mod':
                    gs = [api.Graph.fromSMILES(s, add=False) for s in seeds]
                    rs = [api.Rule.fromGMLString(s, add=False) for s in gml_rules]
                else:
                    gs = [api.smiles(s, add=False) for s in seeds]
                    rs = [api.ruleGMLString(s, add=False) for s in gml_rules]
                # Match the published no-self-mixtures policy.
                strat = api.addSubset(gs) >> api.repeat[args.depth](api.leftPredicate[
                    lambda d: len(set(g.id for g in d.left)) == len(d.left)
                ](rs))
                dg = api.DG(graphDatabase=gs)
                with dg.build() as builder:
                    builder.execute(strat, verbosity=0)
            elapsed = time.perf_counter() - started
            if repeat:
                result['samples'].append(elapsed)
            if args.tool == 'synkit':
                species = [d['smiles'] for _,d in graph.nodes(data=True) if d.get('kind') == 'species']
                if not species:
                    species = list(builder._species_index)
                reactions = []
                for n,d in graph.nodes(data=True):
                    if n in builder._species_index.values():
                        continue
                    left, right = [], []
                    for u,_,e in graph.in_edges(n,data=True):
                        left.extend([graph.nodes[u].get('smiles',graph.nodes[u].get('label'))] * int(e.get('stoich',e.get('weight',1))))
                    for _,v,e in graph.out_edges(n,data=True):
                        right.extend([graph.nodes[v].get('smiles',graph.nodes[v].get('label'))] * int(e.get('stoich',e.get('weight',1))))
                    if left or right:
                        reactions.append((left,right))
                result['native_steps'] = builder.native_report
            else:
                species = [v.graph.smiles for v in dg.vertices]
                reactions = [([v.graph.smiles for v in e.sources], [v.graph.smiles for v in e.targets]) for e in dg.edges]
            snapshot = normalize(species, reactions)
            if repeat and snapshot != result['snapshot']:
                raise AssertionError('Network changed between repeats')
            result['snapshot'] = snapshot
            result['n_species'],result['n_reactions'] = len(snapshot['species']),len(snapshot['reactions'])
            result['status'] = 'running'
            Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
            print(json.dumps({'repeat':repeat,'depth':args.depth,'tool':args.tool,'seconds':elapsed,
                              'species':result['n_species'],'reactions':result['n_reactions']}), flush=True)
    result['status'] = 'complete'
    Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tool', choices=['synkit','mod','molfoundry'],required=True)
    p.add_argument('--depth', type=int,required=True)
    p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--inputs',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--synkit-root',default=os.environ.get('SYNKIT_ROOT',''))
    p.add_argument('--expander')
    run(p.parse_args())
