"""Render and verify the genuine DPO comparison, separately from RDKit controls."""
import hashlib
import json
from pathlib import Path

PAPER=Path(__file__).resolve().parents[2]
EVIDENCE=PAPER/'analysis/formose-expansion'


def read(name):
    return json.loads((EVIDENCE/name).read_text())


def primary_report():
    path=EVIDENCE/'dpo/comparison.json'
    if not path.exists():return None
    report=json.loads(path.read_text())
    return report if any(r['backend']=='synreactor-native' for r in report['runs']) else None


def cell(row):
    if row is None:return '--'
    if row['status']=='complete':
        seconds=row['median_seconds']
        return f"{seconds:.2f}" if seconds>=1 else f"{seconds:.3g}"
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


def macros():
    report=primary_report()
    if report is None:return '% Full native DPO comparison is still running.\n'
    lines=['% Generated from formose-expansion/dpo/comparison.json.']
    for depth,word in [(4,'Four'),(5,'Five'),(6,'Six')]:
        row=next((r for r in report['runs'] if r['depth']==depth and r['backend']=='synreactor-native' and r['status']=='complete'),None)
        if row:lines.append(r'\newcommand{\DpoNative'+word+'Seconds}{'+f"{row['median_seconds']:.3g}"+'}')
    for backend,name in [('mod','Mod'),('molfoundry','Molfoundry')]:
        row=next((r for r in report['runs'] if r['depth']==6 and r['backend']==backend and r['status']=='complete'),None)
        if row:lines.append(r'\newcommand{\Dpo'+name+'SixSeconds}{'+f"{row['median_seconds']:.3g}"+'}')
    return '\n'.join(lines)+'\n'


def verify():
    report=primary_report()
    assert report is not None
    assert report['source_unchanged'] is True
    assert report['repeats_requested']==5 and report['warmups']==1
    original=read('comparison.json')
    checks=[]
    for row in report['runs']:
        if row['status']!='complete':continue
        data=read('dpo/'+row['result'])
        assert len(data['samples'])==report['repeats_requested'] and data['status']=='complete'
        snap=data['snapshot'].copy();digest=snap.pop('sha256')
        assert hashlib.sha256(json.dumps(snap,sort_keys=True).encode()).hexdigest()==digest
        peers=[r for r in original['runs'] if r['depth']==row['depth'] and r['status']=='complete']
        assert all(r['snapshot_sha256']==digest for r in peers)
        identity=data.get('dpo_report')
        if row['backend']=='synreactor-native':
            assert identity['controller']=='cpp' and 'dpo' in identity['engine']
        checks.append({'depth':row['depth'],'network_verified':True,'engine':identity['engine'] if identity else row['backend'],
                       'controller':identity['controller'] if identity else 'external'})
    (EVIDENCE/'dpo/verification.json').write_text(json.dumps({'completed_networks':checks},indent=2)+'\n')
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
                        textcoords='offset points',ha='right',va='center',fontsize=8)
    growth.set(yscale='log',xlabel='Rule-expansion depth',ylabel='Network size',ylim=(.7,25000))
    growth.legend(fontsize=8,loc='upper left')
    for offset,(backend,label,marker,style,color) in zip([-.05,0,.05],[
        ('synreactor-native','SynReactor C++ DPO','o','-','black'),
        ('mod','MØD 1.0.0.75','s','--','#555555'),
        ('molfoundry','MolFoundry 0.1.11','^',':','#888888')]):
        rows=[r for r in report['runs'] if r['backend']==backend and r['status']=='complete']
        for row in rows:
            samples=row['samples']
            ax.scatter([row['depth']+offset+(i-2)*.012 for i in range(5)],samples,
                       marker=marker,s=9,color=color,alpha=.35,zorder=2,linewidths=0)
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
    fig.savefig(EVIDENCE/'dpo/runtime.png',dpi=180)
    plt.close(fig)
    print(f'Verified {len(checks)} completed generation cases, each with five timed repeats.')


if __name__=='__main__':verify()
