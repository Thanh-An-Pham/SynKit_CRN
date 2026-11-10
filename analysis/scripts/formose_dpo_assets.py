"""Render and audit native SynReactor, MØD and MolFoundry formose benchmarks."""
import json
from pathlib import Path
import statistics

PAPER=Path(__file__).resolve().parents[2]
EVIDENCE=PAPER/'analysis/formose-expansion'


def read(name):
    return json.loads((EVIDENCE/name).read_text())


def primary_report():
    path=EVIDENCE/'dpo/comparison.json'
    if not path.exists():return None
    report=json.loads(path.read_text())
    if report.get('schema')=='synkit.formose-dpo-comparison/5':
        required={(depth,backend) for depth in range(1,7)
                  for backend in ['synreactor-native','mod','molfoundry']}
        if not (report.get('status')=='PASS' and report.get('complete') is True
                and report.get('cross_tool_equality_checked') is True
                and report.get('repeats_requested')==20 and report.get('warmups')==3
                and len(report['runs'])==18
                and {(r['depth'],r['backend']) for r in report['runs']}==required
                and all(r['status']=='complete' and len(r['samples'])==20
                        and len(r['warmup_samples'])==3 for r in report['runs'])):
            raise ValueError('The manuscript requires the complete matched six-depth, 20-repeat DPO comparison')
    return report if any(r['backend']=='synreactor-native' for r in report['runs']) else None


def cell(row):
    if row is None:return '--'
    if row['status']=='complete':
        seconds=row['median_seconds']
        return f"{seconds:.3g}"
    if row['status']=='timeout':return 'timeout'
    path=EVIDENCE/row['log'] if row['log'].startswith('runs/') else EVIDENCE/'dpo'/row['log']
    log=path.read_text()
    if 'binding_limit' in log:return 'binding limit'
    return 'memory' if 'MemoryError' in log or 'oom' in log.lower() else 'error'


def table():
    report=primary_report()
    if report is None:return '% Full native DPO comparison is still running.\n'
    lines=[]
    for depth in sorted({r['depth'] for r in report['runs']}):
        rows={r['backend']:r for r in report['runs'] if r['depth']==depth}
        backends=['synreactor-native','mod','molfoundry']
        if not all(name in rows and rows[name]['status']=='complete' for name in backends):
            continue
        known=rows['synreactor-native']
        values=[str(depth),str(known['n_species']),str(known['n_reactions'])]
        values.extend(cell(rows[name]) for name in backends)
        lines.append(' & '.join(values)+r' \\')
    return '\n'.join(lines)+'\n'


def verify(*, render=True):
    report=primary_report()
    assert report is not None
    assert report['source_unchanged'] is True
    assert report['repeats_requested']==20 and report['warmups']==3
    assert report['timeout_seconds_per_expansion']==180
    assert report['schema']=='synkit.formose-dpo-comparison/5'
    assert report['status']=='PASS' and report['complete'] is True
    assert report['cross_tool_equality_checked'] is True
    assert all(report['checks'].values())
    assert len(report['runs'])==18
    assert {(r['depth'],r['backend']) for r in report['runs']} == {
        (depth,backend) for depth in range(1,7)
        for backend in ['synreactor-native','mod','molfoundry']}
    checks=[]
    snapshots={}
    from rdkit import Chem
    for row in report['runs']:
        assert row['status']=='complete', (row['backend'],row['depth'],row['status'])
        data=read('dpo/'+row['result'])
        assert len(data['samples'])==report['repeats_requested'] and data['status']=='complete'
        assert len(data['warmup_samples'])==report['warmups']
        assert len(data['audits'])==report['repeats_requested']+report['warmups']
        assert [(a['phase'],a['index']) for a in data['audits']]==(
            [('warmup',i) for i in range(report['warmups'])]+
            [('measured',i) for i in range(report['repeats_requested'])])
        assert data['samples']==row['samples'] and data['warmup_samples']==row['warmup_samples']
        assert all(value>0 for value in data['samples'])
        q=statistics.quantiles(data['samples'],n=4,method='inclusive')
        assert (row['median_seconds'],row['q1_seconds'],row['q3_seconds'])==(statistics.median(data['samples']),q[0],q[2])
        snap=data['snapshot']
        assert data['snapshot_verified'] is True
        assert all(audit['snapshot_matches'] for audit in data['audits'])
        assert snapshots.setdefault(row['depth'], snap) == snap
        events=[event for event in report['schedule'] if event['depth']==row['depth'] and event['backend']==row['backend']]
        assert [(event['phase'],event['index']) for event in events]==[(audit['phase'],audit['index']) for audit in data['audits']]
        assert all(event['snapshot_matches'] for event in events)
        assert [event['seconds'] for event in events if event['phase']=='measured']==data['samples']
        assert row['snapshot_matches'] is True
        assert (row['n_species'],row['n_reactions'])==(len(snap['species']),len(snap['reactions']))
        peers=[r for r in report['runs'] if r['depth']==row['depth'] and r['status']=='complete']
        assert {'mod','molfoundry'} <= {peer['backend'] for peer in peers}
        assert all(r['snapshot_matches'] for r in peers)
        carbons={s:sum(atom.GetAtomicNum()==6 for atom in Chem.MolFromSmiles(s).GetAtoms()) for s in snap['species']}
        for left,right in snap['reactions']:
            assert sum(carbons[s]*n for s,n in left)==sum(carbons[s]*n for s,n in right)
        identity=data.get('dpo_report')
        if row['backend']=='synreactor-native':
            assert identity['controller']=='cpp' and 'dpo' in identity['engine']
            assert identity['status']=='complete' and identity['stop_reason']=='repeats'
            assert identity['steps_run']==row['depth']
            assert all(not step['capped_rules'] and not step['task_budget_exhausted'] for step in identity['steps'])
        if row['backend']=='synreactor-native':
            assert identity['hydrogen_model']=='implicit' and identity['component_policy']=='union'
            assert identity['radical_policy']=='ignore'
            assert identity['engine']=='synkit_dpo_cpp'
        checks.append({'depth':row['depth'],'backend':row['backend'],'network_verified':True,
                       'carbon_balance_verified':True,
                       'warmups_verified':len(data['warmup_samples']),'measurements_verified':len(data['samples']),
                       'engine':identity['engine'] if identity else row['backend'],
                       'controller':identity['controller'] if identity else 'external'})
    verification={
        'completed_networks':checks,
        'verified_expansions':sum(c['warmups_verified']+c['measurements_verified'] for c in checks),
        'failures':[r for r in report['runs'] if r['status']!='complete'],
    }
    if not render:
        return verification
    (EVIDENCE/'dpo/verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    from result_style import configure, plain_heading
    import matplotlib.pyplot as plt
    configure()
    fig,(growth,ax)=plt.subplots(1,2,figsize=(7.2,3.0))
    native_rows=[r for r in report['runs'] if r['backend']=='synreactor-native' and r['status']=='complete']
    for key,label,marker,face,style,offset in [
        ('n_species','Species','o','black','-',-13),
        ('n_reactions','Reactions','s','white','--',8)]:
        growth.plot([r['depth'] for r in native_rows],[r[key] for r in native_rows],
                    marker=marker,mfc=face,mec='black',color='black',linestyle=style,label=label,ms=4)
        last=native_rows[-1]
        growth.annotate(f"{last[key]:,}",(last['depth'],last[key]),xytext=(-5,offset),
                        textcoords='offset points',ha='right',va='center',fontsize=8,
                        bbox={'facecolor':'white','edgecolor':'none','pad':.2})
    growth.set(yscale='log',xlabel='Rule-expansion depth',ylabel='Network size',ylim=(.7,40000))
    growth.legend(fontsize=8,loc='upper left')
    for offset,(backend,label,marker,style,color) in zip([-.05,0,.05],[
        ('synreactor-native','SynReactor C++ DPO','o','-','black'),
        ('mod','MØD 1.0.0.75','s','--','#555555'),
        ('molfoundry','MolFoundry 0.1.11','^',':','#888888')]):
        rows=[r for r in report['runs'] if r['backend']==backend and r['status']=='complete']
        for row in rows:
            samples=row['samples']
            ax.scatter([row['depth']+offset+(i-(len(samples)-1)/2)*.006 for i in range(len(samples))],samples,
                       marker=marker,s=7,color=color,alpha=.22,zorder=2,linewidths=0)
        ax.errorbar([r['depth'] for r in rows],[r['median_seconds'] for r in rows],
                    yerr=[[r['median_seconds']-r['q1_seconds'] for r in rows],[r['q3_seconds']-r['median_seconds'] for r in rows]],
                    marker=marker,linestyle=style,color=color,label=label,capsize=2)
    depths=sorted({r['depth'] for r in report['runs'] if r['status']=='complete'})
    ax.set(yscale='log',xlabel='Rule-expansion depth',ylabel='Full API runtime from seeds (s)',
           xlim=(min(depths)-.3,max(depths)+.3))
    for panel,letter in [(growth,'a'),(ax,'b')]:
        panel.set_xticks(depths)
        panel.set_xlim(min(depths)-.3,max(depths)+.3)
        plain_heading(panel,letter,y=1.03)
    ax.set_ylabel('Full API runtime (s)')
    ax.legend(fontsize=7.5,loc='upper left',handlelength=1.8)
    fig.subplots_adjust(left=.085,right=.985,bottom=.19,top=.88,wspace=.40)
    fig.savefig(PAPER/'Figure/formose_dpo.pdf')
    plt.close(fig)
    print(f"Verified {len(checks)} completed generation cases, each with {report['repeats_requested']} timed repeats and {report['warmups']} warm-ups.")


if __name__=='__main__':verify()
