"""Audit saved expansion comparisons and render manuscript table/figure."""
import hashlib
import json
from pathlib import Path
import statistics

from result_style import configure, heading
import matplotlib.pyplot as plt


def failure_label(row):
    if row['status']=='timeout': return 'timeout'
    root=Path(__file__).resolve().parents[2]/'analysis/formose-expansion'
    log=(root/row['log']).read_text()
    if 'MemoryError' in log: return 'memory limit'
    if 'oom' in log.lower(): return 'OOM crash'
    return 'error'


def variant_table(original, variants):
    rows=[]
    for depth in sorted({r['depth'] for r in original['runs']}):
        values=[str(depth)]
        for mode in ['python-explicit','native-explicit','native-lewis']:
            if mode=='native-lewis':
                row=next(r for r in original['runs'] if r['depth']==depth and r['tool']=='synkit')
            else:
                row=next((r for r in variants['runs'] if r['depth']==depth and r['mode']==mode),None)
            if row is None:
                values.append('pending')
                continue
            if row['status']=='complete':
                values.append(f"{row['median_seconds']:.4g} [{row['q1_seconds']:.4g}, {row['q3_seconds']:.4g}]")
            else: values.append(failure_label(row))
        rows.append(' & '.join(values)+r' \\')
    return '\n'.join(rows)+'\n'


def control_summary(original, variants):
    rows={r['mode']:r for r in variants['runs'] if r['depth']==6 and r['status']=='complete'}
    if len(rows)!=3: return '% Implementation controls are still running.\n'
    baseline=next(r for r in original['runs'] if r['tool']=='synkit' and r['depth']==6)
    lewis=baseline['median_seconds'];explicit=rows['native-explicit']['median_seconds'];python=rows['python-explicit']['median_seconds']
    instrumented=rows['native-lewis']
    return (f"At depth six, median generation times are {python:.3g}~s for Python/RDKit, {explicit:.3g}~s for RDKit C++ with explicit hydrogen, "
            f"and {lewis:.3g}~s for RDKit C++ with Lewis hydrogen counts. The RDKit C++ Lewis path is {explicit/lewis:.2f} times faster than "
            f"the RDKit C++ explicit-hydrogen path and {python/lewis:.2f} times faster than the Python/RDKit path. "
            f"A separate instrumented RDKit C++ Lewis control has a median total time of {instrumented['median_seconds']:.3g}~s, "
            f"with {instrumented['median_native_request_seconds']:.3g}~s in the native request, including transport and decoding, "
            f"and {instrumented['median_python_construction_seconds']:.3g}~s in Python-side network construction. "
            "This measurement concerns the RDKit SMARTS controller; SynReactor DPO execution is a separate backend.\n")


def table(report):
    rows=[]
    for depth in sorted({r['depth'] for r in report['runs']}):
        runs={r['tool']:r for r in report['runs'] if r['depth']==depth}
        complete=[r for r in runs.values() if r['status']=='complete']
        counts=([str(complete[0]['n_species']),str(complete[0]['n_reactions'])] if complete else ['--','--'])
        values=[str(depth),*counts]
        for tool in ['synkit','mod','molfoundry']:
            r=runs[tool]
            if r['status']=='complete':
                values.append(f"{r['median_seconds']:.4g} [{r['q1_seconds']:.4g}, {r['q3_seconds']:.4g}]")
            else:
                values.append(failure_label(r))
        rows.append(' & '.join(values)+r' \\')
    return '\n'.join(rows)+'\n'


def main():
    paper=Path(__file__).resolve().parents[2]
    out=paper/'analysis/formose-expansion'
    report=json.loads((out/'comparison.json').read_text())
    for name,digest in report['source_sha256'].items():
        source=out/name if name=='inputs.json' else Path(__file__).with_name(name)
        assert hashlib.sha256(source.read_bytes()).hexdigest()==digest,(name,'benchmark source changed')
    checks=[]
    for row in report['runs']:
        if row['status']!='complete':
            continue
        data=json.loads((out/row['result']).read_text())
        assert data['status']=='complete'
        assert len(data['samples'])==report['repeats_requested']
        assert statistics.median(data['samples'])==row['median_seconds']
        snapshot=data['snapshot'].copy();digest=snapshot.pop('sha256')
        assert hashlib.sha256(json.dumps(snapshot,sort_keys=True).encode()).hexdigest()==digest
        assert digest==row['snapshot_sha256']
        assert len(snapshot['species'])==row['n_species']
        assert len(snapshot['reactions'])==row['n_reactions']
        # Check carbon balance independently from the exported incidence multisets.
        from rdkit import Chem
        carbons={s:sum(a.GetAtomicNum()==6 for a in Chem.MolFromSmiles(s).GetAtoms()) for s in snapshot['species']}
        for left,right in snapshot['reactions']:
            assert sum(carbons[s]*n for s,n in left)==sum(carbons[s]*n for s,n in right)
        if row['tool']=='synkit':
            native=data['native_steps']
            assert native['steps_run']==row['depth']
            assert all(not s['capped_rules'] and not s['task_budget_exhausted'] for s in native['steps'])
        checks.append({'tool':row['tool'],'depth':row['depth'],'snapshot_hash_and_carbon_balance_ok':True})
    for depth,eq in report['equivalence'].items():
        complete=[r for r in report['runs'] if r['depth']==int(depth) and r['status']=='complete']
        if len(complete)>1:
            assert eq['identical_snapshots']
    routes=json.loads((out/'route-comparison.json').read_text())
    assert routes['all_ten_agree'] and routes['positive_certificates_replayed']==3
    controls=[]
    variants_path=out/'variants.json'
    variants=json.loads(variants_path.read_text()) if variants_path.exists() else None
    if variants and len(variants['runs'])==15:
        for name,digest in variants['source_sha256'].items():
            assert hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()==digest
        for row in variants['runs']:
            if row['status']!='complete': continue
            data=json.loads((out/row['result']).read_text())
            assert data['status']=='complete' and len(data['samples'])==5
            snapshot=data['snapshot'].copy();digest=snapshot.pop('sha256')
            assert hashlib.sha256(json.dumps(snapshot,sort_keys=True).encode()).hexdigest()==digest
            peer=next(r for r in report['runs'] if r['tool']=='synkit' and r['depth']==row['depth'])
            assert peer['snapshot_sha256']==digest and row['matches_original_native_lewis']
            if row['mode'].startswith('native'):
                assert data['native_report']['steps_run']==row['depth']
                assert all(not s['capped_rules'] and not s['task_budget_exhausted'] for s in data['native_report']['steps'])
            controls.append({'mode':row['mode'],'depth':row['depth'],'matches_verified_native_lewis_snapshot':True})
        (out/'variants-table.tex').write_text(variant_table(report,variants))
    (out/'verification.json').write_text(json.dumps({'complete_network_audits':checks,'variant_audits':controls,'route_comparison_ok':True},indent=2)+'\n')
    (out/'table.tex').write_text(table(report))
    configure()
    fig,axes=plt.subplots(1,2,figsize=(7.2,2.8))
    for tool,marker,style,color,label in [('synkit','o','-','black','SynKit RDKit C++ (Lewis H)'),('mod','s','--','#555555','MØD 1.0.0.75'),('molfoundry','^',':','#888888','MolFoundry 0.1.11')]:
        rows=[r for r in report['runs'] if r['tool']==tool and r['status']=='complete']
        axes[0].errorbar([r['depth'] for r in rows],[r['median_seconds'] for r in rows],
                        yerr=[[r['median_seconds']-r['q1_seconds'] for r in rows],[r['q3_seconds']-r['median_seconds'] for r in rows]],
                        marker=marker,linestyle=style,color=color,label=label,capsize=2)
    axes[0].set_yscale('log');axes[0].set_ylabel('Cumulative generation time (s)');axes[0].legend(fontsize=7)
    rows=[r for r in report['runs'] if r['tool']=='synkit' and r['status']=='complete']
    if variants and len(variants['runs'])==15:
        for mode,marker,style,color,label in [('python-explicit','s','--','#888888','Python/RDKit, explicit H'),('native-explicit','^',':','#555555','RDKit C++, explicit H'),('native-lewis','o','-','black','RDKit C++, Lewis H')]:
            selected=rows if mode=='native-lewis' else [r for r in variants['runs'] if r['mode']==mode and r['status']=='complete']
            axes[1].errorbar([r['depth'] for r in selected],[r['median_seconds'] for r in selected],
                            yerr=[[r['median_seconds']-r['q1_seconds'] for r in selected],[r['q3_seconds']-r['median_seconds'] for r in selected]],
                            marker=marker,linestyle=style,color=color,label=label,capsize=2)
        axes[1].set_ylabel('Cumulative generation time (s)')
    else:
        for key,marker,style,label in [('n_species','o','-','Species'),('n_reactions','s','--','Reactions')]:
            axes[1].plot([r['depth'] for r in rows],[r[key] for r in rows],marker=marker,linestyle=style,color='black',label=label)
        axes[1].set_ylabel('Network size')
    axes[1].set_yscale('log');axes[1].legend(fontsize=7)
    for ax in axes:
        ax.set_xlabel('Rule-expansion depth');ax.set_xticks(range(1,8));ax.set_xlim(.7,7.3)
    heading(axes[0],'a','Expansion runtime');heading(axes[1],'b','SynKit implementation and H model' if variants and len(variants['runs'])==15 else 'Network growth')
    fig.subplots_adjust(left=.085,right=.98,bottom=.19,top=.84,wspace=.35)
    fig.savefig(paper/'Figure/formose_expansion.pdf')
    plt.close(fig)
    print(f'Verified {len(checks)} completed networks and ten route decisions; generated table and figure.')


if __name__=='__main__':
    main()
