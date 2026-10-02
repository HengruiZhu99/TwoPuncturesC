"""Four system/Krylov combinations plus archived default-state controls.

Fresh workers, one numerical job at a time. Computational norms are matched;
physical accuracy is compared only within each equation system.
"""
import argparse,json,os,platform,re,subprocess,sys
from pathlib import Path
import numpy as np
from benchmark_bowen_york import ROOT,digest
from benchmark_common_stopping import compare


def arrays(a,b):
    return {k:dict(finite=bool(np.isfinite(b[k]).all()),bitwise_equal=a[k].tobytes()==b[k].tobytes(),
                   max_scaled_difference=float(np.max(np.abs(a[k]-b[k])/(1+np.abs(a[k]))))) for k in a.files}


def by_compare(reference,candidate):
    result=compare(reference,candidate)
    result['all_retained_arrays_finite']=all(v['finite'] for v in result['arrays'].values())
    result['passed']=result['passed'] and result['all_retained_arrays_finite']
    return result


def hi_coefficients(values,grid):
    na,nb,nphi=grid;v=values.reshape(nphi,nb,na,4)
    def matrix(n):
        i=np.arange(n)[:,None];j=np.arange(n)[None,:]
        return (2/n)*(-1.)**i*np.cos(np.pi*i*(j+.5)/n)
    return np.einsum('ai,bj,kjiv->kbav',matrix(na),matrix(nb),v,optimize=True)


def hi_compare(reference,candidate):
    a=np.load(ROOT/reference['state']);b=np.load(ROOT/candidate['state']);metrics=arrays(a,b)
    ac=hi_coefficients(a['unknowns'],reference['config']['n']);bc=hi_coefficients(b['unknowns'],candidate['config']['n'])
    coefficients=float(np.max(np.abs(ac-bc)/(1+np.abs(ac))))
    core=('unknowns','lapse','gamma','Kij','psi','correction','conformal_metric','Atilde','mean_curvature','dgamma')
    fields=all(metrics[k]['finite'] and metrics[k]['max_scaled_difference']<=1e-10 for k in core)
    charges=float(np.max(np.abs(np.asarray(reference['charges'])-candidate['charges'])/(1+np.abs(reference['charges']))))
    physical={}
    for k in ('H','M'):
        ar=[np.asarray(r[k]) for r in reference['physical_check']['records']]
        br=[np.asarray(r[k]) for r in candidate['physical_check']['records']]
        floor=max(float(np.max(np.abs(x-y))) for series in (ar,br) for x,y in zip(series,series[1:]))
        delta=max(float(np.max(np.abs(x-y))) for x,y in zip(ar,br));bound=max(1e-7,5*floor)
        physical[k]=dict(max_difference=delta,bound=bound,passed=delta<=bound)
    residual=max(float(np.max(np.abs(s['residual']))) for s in (a,b))
    passed=all(v['finite'] for v in metrics.values()) and fields and coefficients<=1e-10 and charges<=1e-10 and residual<=1e-12 and all(p['passed'] for p in physical.values())
    return dict(passed=bool(passed),arrays=metrics,coefficients_scaled_difference=coefficients,charges_scaled_difference=charges,physical=physical,max_recomputed_weighted_F=residual)


def checks(mode,r):
    n=r['diagnostics']['newton_iterations'] if mode=='hispid' else r['work_statistics']['newton_iterations']
    history=r['linear_history'];c=dict(converged=r['converged'],scaling=r['residual_scaling']=='sin3_alpha_beta',
          residual=max(r['computational_norms']['native_weighted']['linf'])<=1e-12,history_complete=len(history)==n)
    if mode=='hispid':
        c['linear']=all(np.isfinite(row).all() and row[1]==.001 and row[2]<=.001 for row in map(np.asarray,history))
        c['iterations']=sum(int(row[3]) for row in history)==r['diagnostics']['krylov_iterations']
    else:
        c['linear']=all(np.isfinite(list(row.values())).all() and row['true_relative_l2']<=.001 and row['true_l2']<=row['absolute_target'] for row in history)
        c['failures']=r['work_statistics']['linear_failures']==r['work_statistics']['modal_failures']==0
        c['iterations']=sum(r['linear_iterations'])==r['work_statistics']['krylov_iterations']
    c['passed']=all(c.values());return c


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('hi6','hi3','by-library','output'):p.add_argument('--'+name,required=True)
    p.add_argument('--repeats',type=int,default=2);args=p.parse_args()
    out=Path(args.output);raw=ROOT/'validation/raw'/out.stem
    if out.exists():raise FileExistsError(out)
    if args.repeats<2:p.error('two interleaved repeats required')
    raw.mkdir(parents=True,exist_ok=False)
    old=json.loads((ROOT/'validation/common_stopping_moderate_40.json').read_text())
    input_path=raw/'input.json';input_path.write_text(json.dumps(old['records']['hi_native'][0]['config'],indent=2)+'\n')
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    def run(label,mode,krylov=None):
        result=raw/(label+'.json');state=raw/(label+'.npz');log=raw/(label+'.log')
        cmd=[sys.executable,str(ROOT/'validation/benchmark_bowen_york.py'),'--hispid-library',args.hi3 if krylov else args.hi6,
             '--by-library',args.by_library,'--case','moderate','--grid','40:80:16','--tolerance','1e-12','--output',str(out),
             '--mode',mode,'--worker-output',str(result),'--state-output',str(state)]
        if mode=='by':cmd+=['--input',str(input_path),'--by-verbose']
        if krylov:
            cmd+=['--krylov',krylov,'--linear-rtol','.001','--physical-check']
            if mode=='by':cmd+=['--by-preconditioner','1','--krylov-maxit','2000']
        print('Starting',label,flush=True)
        with log.open('w') as stream:subprocess.run(cmd,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=1200)
        r=json.loads(result.read_text());r.update(state=str(state.relative_to(ROOT)),state_sha256=digest(state),log=str(log.relative_to(ROOT)),log_sha256=digest(log),command=cmd)
        if mode=='by':
            text=log.read_text()
            r['linear_history']=[dict(true_l2=float(v),true_relative_l2=float(rel),absolute_target=float(target)) for v,rel,target in re.findall(r'linear_true: ([^\s]+) relative ([^\s]+) target ([^\s]+)',text)]
            r['linear_iterations']=[max(map(int,re.findall(r'^(?:bicgstab|gmres):\s+(\d+)\s',part,re.M)),default=0) for part in re.split(r'(?:bicgstab|gmres):  itmax.*\n',text)[1:]]
        print(label,round(r['solve_seconds'],6),r['max_rss_bytes'],r['converged'],flush=True)
        return r
    defaults={}
    for mode,label in (('hispid','hi_native'),('by','by_native')):
        reference=old['records'][label][0];r=run(label+'_default',mode)
        if digest(ROOT/reference['state'])!=reference['state_sha256']:raise RuntimeError('archived reference hash changed')
        a=np.load(ROOT/reference['state']);b=np.load(ROOT/r['state']);metrics=arrays(a,b)
        counts=('newton_iterations','krylov_iterations') if mode=='hispid' else ('newton_iterations','krylov_iterations','jvp_applications','preconditioner_applications','relaxation_sweeps')
        refstats=reference['diagnostics'] if mode=='hispid' else reference['work_statistics']
        stats=r['diagnostics'] if mode=='hispid' else r['work_statistics']
        equal_counts=all(stats[k]==refstats[k] for k in counts)
        trace=lambda record:[line for line in (ROOT/record['log']).read_text().splitlines() if line.startswith(('Newton:','bicgstab:'))]
        trace_equal=mode=='hispid' or trace(r)==trace(reference)
        passed=r['converged'] and set(a.files)==set(b.files) and all(v['bitwise_equal'] and v['finite'] for v in metrics.values()) and equal_counts and trace_equal
        defaults[mode]=dict(passed=bool(passed),arrays=metrics,counts_equal=equal_counts,trace_equal=trace_equal,record=r,reference_state_sha256=reference['state_sha256'])
    protocols=[('hi_gmres','hispid','gmres'),('hi_bicgstab','hispid','bicgstab'),('by_gmres','by','gmres'),('by_bicgstab','by','bicgstab')]
    records={label:[] for label,_,_ in protocols}
    for repeat in range(args.repeats):
        for label,mode,krylov in protocols if repeat%2==0 else reversed(protocols):
            r=run(f'{label}_{repeat}',mode,krylov);r['checks']=checks(mode,r);records[label].append(r)
    equivalence={}
    for prefix,comparison in (('hi',hi_compare),('by',by_compare)):
        ref=records[prefix+'_gmres'][0]
        equivalence[prefix]={label:[comparison(ref,r) for r in series] for label,series in records.items() if label.startswith(prefix+'_')}
    performance={label:{key:float(np.median([r[key] for r in series])) for key in ('solve_seconds','ready_to_sample_seconds','max_rss_bytes')} for label,series in records.items()}
    passed=all(d['passed'] for d in defaults.values()) and all(r['checks']['passed'] for series in records.values() for r in series) and all(c['passed'] for group in equivalence.values() for series in group.values() for c in series)
    result=dict(passed=bool(passed),default_controls=defaults,records=records,equivalence=equivalence,performance=performance,
          repeats=args.repeats,grid=[40,80,16],threads=1,criteria_sha256=digest(ROOT/'validation/krylov_matrix_acceptance.json'),platform=platform.platform(),python=sys.version,
          source_sha256={str(f.relative_to(ROOT)):digest(f) for d in ('src','include','python') for f in (ROOT/d).glob('*') if f.is_file()},
          note='True RHS-relative L2 .001; cubic outer max1e-12; modal M fixed per system; maxK2000/restart64. Only within-system data equivalence; no high-spin binary acceptance inherited.')
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(performance,indent=2),flush=True)
    if not passed:raise SystemExit(1)
if __name__=='__main__':main()
