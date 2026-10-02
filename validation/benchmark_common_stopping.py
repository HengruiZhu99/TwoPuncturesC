"""Matched computational stopping norms, sequential fresh workers and BY data controls.

Only one worker runs at a time. Physical maps/seed geometries are not equated.
"""
import argparse,json,os,platform,re,subprocess,sys
from pathlib import Path
import numpy as np
from benchmark_bowen_york import ROOT,digest


def compare(reference,candidate):
    a=np.load(ROOT/reference['state']);b=np.load(ROOT/candidate['state'])
    metrics={k:dict(finite=bool(np.isfinite(b[k]).all()),bitwise_equal=a[k].tobytes()==b[k].tobytes(),
                   max_absolute_difference=float(np.max(np.abs(a[k]-b[k]))),
                   max_scaled_difference=float(np.max(np.abs(a[k]-b[k])/(1+np.abs(a[k]))))) for k in a.files}
    core=('v_d0','u_d0','cf_v_d0','spectral_v_samples','lapse','psi','gamma','K')
    fields=all(metrics[k]['finite'] and metrics[k]['max_scaled_difference']<=1e-10 for k in core)
    charges={k:float(np.max(np.abs(np.asarray(reference['diagnostics'][k])-candidate['diagnostics'][k])/(1+np.abs(reference['diagnostics'][k])))) for k in ('adm_energy','puncture_masses')}
    physical={}
    for k in ('H','M'):
        ar=[np.asarray(r[k]) for r in reference['physical_check']['records']]
        br=[np.asarray(r[k]) for r in candidate['physical_check']['records']]
        floor=max(float(np.max(np.abs(x-y))) for series in (ar,br) for x,y in zip(series,series[1:]))
        delta=max(float(np.max(np.abs(x-y))) for x,y in zip(ar,br));bound=max(1e-7,5*floor)
        physical[k]=dict(max_difference=delta,step_variation=floor,bound=bound,passed=delta<=bound)
    recomputed=max(float(np.max(np.abs(s['recomputed_F']))) for s in (a,b))
    passed=reference['converged'] and candidate['converged'] and fields and recomputed<=1e-12 and all(v<=1e-10 for v in charges.values()) and all(v['passed'] for v in physical.values())
    return dict(passed=bool(passed),arrays=metrics,charges=charges,physical=physical,max_recomputed_weighted_F=recomputed)


def protocol_checks(label,record):
    hi=label.startswith('hi_');common='common' in label
    stats=record['diagnostics'] if hi else record['work_statistics']
    history=record['linear_history'];n=stats['newton_iterations']
    expected_tag='sin6_alpha_beta' if label=='hi_native' else 'sin3_alpha_beta'
    checks=dict(converged=bool(record['converged']),weight_tag=record['residual_scaling']==expected_tag,
                forcing_option=record['linear_rtol']==(.001 if common else None),
                native_max=bool(np.isfinite(record['computational_norms']['native_weighted']['linf']).all() and max(record['computational_norms']['native_weighted']['linf'])<=1e-12))
    if hi:
        checks['history_complete']=len(history)==n and all(int(row[0])==q for q,row in enumerate(history))
        checks['krylov_total']=sum(int(row[3]) for row in history)==stats['krylov_iterations']
        checks['true_linear_residuals']=all(np.isfinite(row).all() and row[2]<=row[1]*(1+1e-12) and (not common or row[1]==.001) for row in map(np.asarray,history))
    else:
        measured=common or 'modal' in label
        checks['history_complete']=len(record['newton_trace'])==n+1 and (not measured or len(history)==n)
        checks['linear_failures_zero']=stats['linear_failures']==stats['modal_failures']==0
        checks['options']=record['by_integer_parameters']['TP_preconditioner']==int('modal' in label) and record['by_integer_parameters']['TP_linear_relative']==int(common)
        if measured:checks['true_linear_residuals']=all(np.isfinite(list(row.values())).all() and row['true_l2']<=row['absolute_target'] and (not common or row['true_relative_l2']<=.001*(1+1e-12)) for row in history)
    checks['passed']=all(checks.values());return checks


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('hi6','hi3','by-library'):p.add_argument('--'+name,required=True)
    p.add_argument('--output',required=True);p.add_argument('--repeats',type=int,default=2)
    args=p.parse_args();out=Path(args.output);raw=ROOT/'validation/raw'/out.stem
    if args.repeats<2:p.error('at least two interleaved repeats required')
    if out.exists():raise FileExistsError(out)
    raw.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    protocols={
        'hi_native':('hispid',args.hi6,0,None),
        'hi_cubic_adaptive':('hispid',args.hi3,0,None),
        'hi_common':('hispid',args.hi3,0,.001),
        'by_native':('by',args.hi6,0,None),
        'by_common_lines':('by',args.hi6,0,.001),
        'by_modal_native':('by',args.hi6,1,None),
        'by_common_modal':('by',args.hi6,1,.001)}
    records={k:[] for k in protocols};input_path=raw/'input.json'
    # Bootstrap config using a first Hi worker, then reverse order on second pass.
    for repeat in range(args.repeats):
        labels=list(protocols) if repeat%2==0 else list(reversed(protocols))
        for label in labels:
            mode,image,precond,forcing=protocols[label];name=f'{label}_{repeat}'
            output=raw/(name+'.json');state=raw/(name+'.npz');log=raw/(name+'.log')
            cmd=[sys.executable,str(ROOT/'validation/benchmark_bowen_york.py'),'--hispid-library',image,
                 '--by-library',args.by_library,'--case','moderate','--grid','40:80:16','--tolerance','1e-12',
                 '--output',str(out),'--mode',mode,'--worker-output',str(output),'--state-output',str(state),'--physical-check']
            if mode=='by':cmd+=['--input',str(input_path),'--by-verbose','--by-preconditioner',str(precond)]
            if forcing is not None:cmd+=['--linear-rtol',str(forcing)]
            print('Starting',name,flush=True)
            with log.open('w') as stream:subprocess.run(cmd,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=1200,check=True)
            r=json.loads(output.read_text());r.update(state=str(state.relative_to(ROOT)),state_sha256=digest(state),
                                                     log=str(log.relative_to(ROOT)),log_sha256=digest(log),command=cmd)
            if mode=='hispid' and not input_path.exists():input_path.write_text(json.dumps(r['config'],indent=2)+'\n')
            if mode=='by':
                text=log.read_text()
                r['newton_trace']=[dict(iteration=int(i),weighted_linf=float(v)) for i,v in re.findall(r'Newton: it=(\d+)\s+\|F\|=([^\s]+)',text)]
                r['linear_history']=[dict(true_l2=float(v),true_relative_l2=float(rel),absolute_target=float(target)) for v,rel,target in re.findall(r'linear_true: ([^\s]+) relative ([^\s]+) target ([^\s]+)',text)]
                r['linear_iterations']=[max(map(int,re.findall(r'^bicgstab:\s+(\d+)\s',part,re.M)),default=0) for part in re.split(r'bicgstab:  itmax.*\n',text)[1:]]
            r['protocol_checks']=protocol_checks(label,r)
            records[label].append(r)
            counts=r['work_statistics'] if mode=='by' else r['diagnostics']
            print(name,round(r['solve_seconds'],6),r['max_rss_bytes'],counts['newton_iterations'],counts['krylov_iterations'],r['converged'],flush=True)
    controls={}
    reference=records['by_native'][0]
    for label in ('by_native','by_common_lines','by_modal_native','by_common_modal'):
        controls[label]=[compare(reference,r) for r in records[label]]
    performance={k:{metric:float(np.median([r[metric] for r in series])) for metric in ('solve_seconds','ready_to_sample_seconds','max_rss_bytes')} for k,series in records.items()}
    passed=all(r['protocol_checks']['passed'] for series in records.values() for r in series) and all(c['passed'] for cs in controls.values() for c in cs)
    summary=dict(records=records,performance=performance,by_data_controls=controls,passed=bool(passed),
                 grid=[40,80,16],tolerance=1e-12,common_linear_relative_l2=.001,repeats=args.repeats,cpu_threads=1,
                 criteria_sha256=digest(ROOT/'validation/by_modal_acceptance.json'),worker_sha256=digest(ROOT/'validation/benchmark_bowen_york.py'),
                 platform=platform.platform(),python=sys.version,
                 source_sha256={str(f.relative_to(ROOT)):digest(f) for d in ('src','include','python') for f in (ROOT/d).glob('*') if f.is_file()},
                 note='Cubic row weighting and fixed RHS-relative L2 forcing define a shared computational convention, not equal physical accuracy. BY solves one Hamiltonian equation, Hi four coupled equations, on different maps/seed geometry. Independent Cartesian constraints are reported separately; no physical validation inherited.')
    out.write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(performance,indent=2),flush=True)
    if not passed:raise SystemExit(1)
if __name__=='__main__':main()
