"""Reproducible one-thread physical validation; fails closed on acceptance.

Pass --stage seeds|moderate|highspin|highboost. Higher stages require the saved
preceding gates; all failures and raw sampled data are retained.
"""
from __future__ import annotations
import argparse,copy,json,os,sys,time,hashlib
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from hispid import Backend,Hole
from physical import constraints,norms,charges,extrapolate
from configs import moderate,hs99uu,highboost,as_dict
from prolong import for_backend

ROOT=Path(__file__).resolve().parents[1]
SOLVER_CONTROLS=('n','memory_limit_mib','tolerance','max_newton','max_krylov','krylov_restart')
def free_data(config):
    return {k:v for k,v in as_dict(config).items() if k not in SOLVER_CONTROLS}

REPORT=ROOT/'validation/results.json'
RAW=ROOT/'validation/raw'

def points(config,horizon_scaled=False):
    # Fixed noncollocation directions/radii, repeated at every resolution.
    directions=np.array([[.73,.31,.61],[-.41,.82,.39],[.22,-.51,.83],[-.69,-.44,.57],[.39,.73,-.56],[.81,-.38,-.45]])
    directions/=np.linalg.norm(directions,axis=1)[:,None]
    near=[];inside=[]
    for h in config.hole:
        if h.mass<=0:continue
        center=np.array(h.center)
        if horizon_scaled:
            chi=np.linalg.norm(list(h.spin))/h.mass**2
            velocity=np.array(h.velocity);v2=velocity@velocity;G=1/np.sqrt(1-v2)
            rh=.5*h.mass*np.sqrt(1-chi**2)
            ray_horizon=rh/np.sqrt(1+(G*G-1)*(directions@(velocity/np.sqrt(v2) if v2 else np.zeros(3)))**2)
            for factor in (1.5,3):near.extend(center+factor*ray_horizon[:,None]*directions)
            radii=(.6*h.mass,h.mass)
        else:radii=(.6*h.mass,h.mass,1.6*h.mass)
        for radius in radii:near.extend(center+radius*directions)
    for h in range(2):
        rmax=config.inner_max[h]
        if rmax>0:
            for radius in (.3,.7,.9):inside.extend(np.array(config.hole[h].center)+radius*rmax*directions)
    bulk=np.r_[directions*4,directions*8,directions*16]
    return np.r_[np.array(near),bulk,np.array(inside).reshape(-1,3)],len(near),len(bulk)

def seeds(backend):
    records=[]
    cases=[('Schwarzschild',Hole(1)),('Kerr',Hole(1,spin=(.2,-.3,.4))),
           ('boosted_Schwarzschild',Hole(1,velocity=(.2,.3,-.1))),
           ('boosted_Kerr',Hole(1,spin=(.2,-.3,.4),velocity=(.2,.3,-.1)))]
    x=np.array([[1.1,.2,.3],[.7,-.9,.4],[2.7,.5,-1.2],[.08,.12,.07]])
    for label,h in cases:
        sample=lambda x:backend.seed(h,x)
        seq=[];start=time.monotonic()
        for step in (.002,.001,.0005):
            r=constraints(sample,x,step);seq.append(dict(step=step,norms=norms(r)))
        radii=[40.,80.,160.,320.]
        q=[charges(sample,R,ntheta=16,nphi=32).tolist() for R in radii];fit=extrapolate(radii,q)
        vel=np.array(h.velocity);spin=np.array(h.spin);G=1/np.sqrt(1-vel@vel)
        expectedJ=G*spin-(G*G/(G+1))*(vel@spin)*vel
        expected=np.r_[G*h.mass,G*h.mass*vel,expectedJ];err=np.abs(fit-expected)
        passed=seq[-1]['norms']['H_rms']<1e-7 and seq[-1]['norms']['M_rms']<1e-7 and np.max(err)<2e-5
        records.append(dict(case=label,hole=as_dict(h),library_sha256=backend.library_sha256(),step_sequence=seq,radii=radii,charges=q,
                            extrapolated=fit.tolist(),expected=expected.tolist(),charge_error=err.tolist(),
                            seconds=time.monotonic()-start,passed=bool(passed)))
        print(label,seq[-1]['norms'],'charge error',err,flush=True)
    return {'records':records,'passed':all(r['passed'] for r in records)}

def solve_case(backend,factory,levels,label,horizon_scaled=False,adaptive_steps=False,previous_records=None,initial_record=None):
    records=list(previous_records or [])
    previous_values=None;previous_shape=None;previous_config=None
    if records or initial_record:
        last=records[-1] if records else initial_record
        if last['library_sha256']!=backend.library_sha256() or last.get('unknown_parameterization_id')!=backend.parameterization():
            raise ValueError('resume requires the same native library SHA')
        previous_shape=last['resolution'];previous_config={k:v for k,v in last['config'].items() if k not in SOLVER_CONTROLS}
        n,_,np_=previous_shape
        previous_values=np.load(RAW/f"{last['case']}_{n}_{np_}.npz")['unknowns']
    for n,nphi in levels:
        cfg=factory(backend,n,nphi);x,near,bulk=points(cfg,horizon_scaled);start=time.monotonic()
        step=np.minimum(.002,.001*np.min([np.linalg.norm(x-np.array(h.center),axis=1) for h in cfg.hole if h.mass>0],axis=0)) if adaptive_steps else np.full(len(x),.002)
        rec=dict(case=label,config=as_dict(cfg),resolution=list(cfg.n),verifier_order=4,verifier_step=.002,
                 library_sha256=backend.library_sha256(),residual_scaling=backend.residual_scaling(),
                 near_sample_count=near,bulk_sample_count=bulk,attenuation_sample_count=len(x)-near-bulk,
                 horizon_scaled=horizon_scaled,verifier_steps=step.tolist(),unknown_parameterization=backend.parameterization_description(),unknown_parameterization_id=backend.parameterization(),collocation_maps=backend.parameterization_maps())
        with backend.create(cfg) as s:
            comparable=free_data(cfg)
            if previous_values is not None and comparable==previous_config:
                s.set_unknowns(for_backend(backend,previous_values,previous_shape,list(cfg.n)))
                rec['initial_guess_from_resolution']=previous_shape
            rec['creation_seconds']=time.monotonic()-start;rec['diagnostics']=s.solve()
            collocation=s.equation_samples();equivalent=collocation['physical_equivalent'];cg=collocation['attenuation']
            rec['physical_equivalent_g1_linf']=np.max(np.abs(equivalent[cg==1]),axis=0).tolist()
            rec['equation_normalized_g_lt_one_linf']=np.max(np.abs(equivalent[cg<1]),axis=0).tolist() if np.any(cg<1) else []
            r=constraints(s.sample,x,step)
            rec['near']=norms(r,np.arange(len(x))<near)
            rec['bulk']=norms(r,(np.arange(len(x))>=near)&(np.arange(len(x))<near+bulk))
            rec['attenuation']=norms(r,r['attenuation']<1)
            rec['g_equals_one']=norms(r,r['attenuation']==1)
            rec['min_metric_eigenvalue']=float(np.min(r['min_metric_eigenvalue']))
            # Record verifier-step sensitivity on near/exterior only.
            rec['verifier_step_check']=[]
            for step in (.004,.001):
                factor=step/.002
                rr=constraints(s.sample,x[:near+bulk],factor*np.array(rec['verifier_steps'][:near+bulk]))
                rec['verifier_step_check'].append(dict(step_factor=factor,near=norms(rr,np.arange(near+bulk)<near),bulk=norms(rr,np.arange(near+bulk)>=near)))
            radii=[100.,200.,400.]
            rec['charge_radii']=radii;rec['charges']=[s.charges(R,ntheta=12,nphi=24).tolist() for R in radii]
            rec['charges_extrapolated']=extrapolate(radii,rec['charges']).tolist()
            RAW.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(RAW/f'{label}_{n}_{nphi}.npz',unknowns=s.unknowns(),points=x,verifier_steps=np.array(rec['verifier_steps']),**r)
            np.savez_compressed(RAW/f'{label}_{n}_{nphi}_collocation.npz',**collocation)
            rec['passed_local']=rec['diagnostics']['status']==0 and rec['min_metric_eigenvalue']>0 and all(rec[k][q]<1e-4 for k in ('near','bulk') for q in ('H_rms','M_rms'))
            rec['passed_strict']=rec['passed_local'] and all(rec[k][q]<1e-6 for k in ('near','bulk') for q in ('H_rms','M_rms')) and all(rec[k][q]<1e-4 for k in ('near','bulk') for q in ('H_max','M_max'))
            previous_values=s.unknowns();previous_shape=list(cfg.n);previous_config=comparable
        rec['total_seconds']=time.monotonic()-start;records.append(rec)
        print(label,rec['resolution'],rec['diagnostics'],rec['near'],rec['bulk'],rec['charges_extrapolated'],flush=True)
        # Preserve progress even if an expensive later run is interrupted.
        previous=json.loads(REPORT.read_text()) if REPORT.exists() else {};previous[label]={'records':records,'passed':False}
        REPORT.write_text(json.dumps(previous,indent=2)+'\n')
    best=records[-1];fine=records[-3:]
    converges=len(fine)>=3 and all(fine[i+1][k][q]<fine[i][k][q] for i in range(2) for k in ('near','bulk') for q in ('H_rms','M_rms'))
    charge_stable=len(records)>=2 and max(abs(np.array(records[-1]['charges_extrapolated'])-np.array(records[-2]['charges_extrapolated'])))<.005
    return dict(records=records,acceptance_resolutions=[r['resolution'] for r in fine],passed=bool(best['passed_local'] and converges and charge_stable),passed_strict=bool(best['passed_strict'] and converges and charge_stable),converges=bool(converges),charge_stable=bool(charge_stable),horizon_enclosure_verified=False)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--library',required=True);parser.add_argument('--stage',choices=['seeds','moderate','highspin','highboost'],required=True)
    parser.add_argument('--levels',default=None,help='e.g. 12:8,20:12,28:16')
    parser.add_argument('--far-radius',type=float,default=None);parser.add_argument('--label',default=None)
    parser.add_argument('--resume',action='store_true',help='append grids to the same case and warm-start from its last compatible checkpoint')
    parser.add_argument('--evaluate-only',action='store_true',help='with --resume, recompute summary gates without another solve')
    parser.add_argument('--initial-from',help='warm-start a separate diagnostic case from the last compatible saved grid of this case')
    parser.add_argument('--initial-compatibility-proof',help='bitwise field/operator witness for an API-only initial-guess migration; fresh solves and gates are still required')
    parser.add_argument('--memory-mib',type=int,default=None,help='explicit per-context allocation budget (default2048,max8192)')
    parser.add_argument('--krylov-restart',type=int);parser.add_argument('--max-krylov',type=int);parser.add_argument('--max-newton',type=int)
    parser.add_argument('--inner-flatten',type=int,choices=[0,1],help='separately labeled interior-operator experiment')
    parser.add_argument('--inner-throat-window',help='lo:hi fractions of the minimum contracted isolated Kerr throat; a screen, not horizon enclosure')
    parser.add_argument('--angular-ratio',type=float,default=1.,help='N_polar=ceil(ratio*N_radial); applies to every saved level')
    args=parser.parse_args()
    if not np.isfinite(args.angular_ratio) or not 1<=args.angular_ratio<=4:
        parser.error('--angular-ratio must be finite and in[1,4]')
    if args.evaluate_only and (not args.resume or args.stage=='seeds'):
        parser.error('--evaluate-only requires a binary stage and --resume')
    if args.initial_compatibility_proof and not args.initial_from:
        parser.error('--initial-compatibility-proof requires --initial-from')
    backend=Backend(args.library);report=json.loads(REPORT.read_text()) if REPORT.exists() else {}
    sha=backend.library_sha256()
    def gate(label):
        result=report.get(label,{})
        return bool(result.get('passed') and result.get('records') and all(r.get('library_sha256')==sha for r in result['records']))
    if args.stage!='seeds' and not gate('seeds'):raise SystemExit('single-seed gate for this library has not passed')
    if args.stage in ('highspin','highboost'):
        moderate_labels=set(('moderate','moderate_far0','moderate_far0_stable'))|{k for k,v in report.items() if isinstance(v,dict) and v.get('stage')=='moderate'}
        if not any(gate(k) for k in moderate_labels):raise SystemExit('moderate convergence gate for this library has not passed')
        if not gate('covariance') or not gate(report['covariance'].get('source_case','')):raise SystemExit('solved coordinate covariance gate for this library has not passed')
    if args.stage=='seeds':result=seeds(backend);label='seeds'
    else:
        label=args.label or args.stage;selected={'moderate':moderate,'highspin':hs99uu,'highboost':highboost}[args.stage]
        def factory(backend,n,nphi):
            config=selected(backend,n,nphi)
            config.n[1]=int(np.ceil(args.angular_ratio*n))
            if args.far_radius is not None:config.far_radius=args.far_radius
            if args.memory_mib is not None:config.memory_limit_mib=args.memory_mib
            if args.inner_flatten is not None:config.inner_flatten=args.inner_flatten
            if args.inner_throat_window:
                lo,hi=map(float,args.inner_throat_window.split(':'))
                if not 0<=lo<hi<1:raise ValueError('throat window must satisfy0<=lo<hi<1')
                for h,hole in enumerate(config.hole):
                    if hole.mass<=0:continue
                    chi=np.linalg.norm(list(hole.spin))/hole.mass**2
                    v2=np.dot(list(hole.velocity),list(hole.velocity))
                    radius=.5*hole.mass*np.sqrt(1-chi**2)*np.sqrt(1-v2)
                    config.inner_min[h]=lo*radius;config.inner_max[h]=hi*radius
            for field in ('krylov_restart','max_krylov','max_newton'):
                if getattr(args,field) is not None:setattr(config,field,getattr(args,field))
            return config
        levels=args.levels or ('24:12,40:20,56:28' if args.stage=='moderate' else '24:8,40:12,56:16')
        levels=[] if args.evaluate_only else [tuple(map(int,v.split(':'))) for v in levels.split(',')]
        previous=report.get(label,{}).get('records') if args.resume else None
        if args.evaluate_only and not previous:raise ValueError('no saved records to evaluate')
        initial=report[args.initial_from]['records'][-1] if args.initial_from else None
        if initial is not None and args.initial_compatibility_proof:
            from checkpoints import restore_equivalent
            restore_equivalent(backend,initial,json.loads(Path(args.initial_compatibility_proof).read_text()))
            initial=copy.deepcopy(initial)
            initial['initial_guess_source_library_sha256']=initial['library_sha256']
            initial['library_sha256']=sha
        result=solve_case(backend,factory,levels,label,horizon_scaled=args.stage in ('highspin','highboost'),adaptive_steps=args.stage in ('highspin','highboost'),previous_records=previous,initial_record=initial)
        if args.initial_compatibility_proof:
            result.update(initial_guess_source_case=args.initial_from,
                initial_guess_source_library_sha256=initial['initial_guess_source_library_sha256'],
                initial_guess_compatibility_proof=args.initial_compatibility_proof)
        if args.stage=='highspin':
            error=abs(result['records'][-1]['charges_extrapolated'][0]-.980124)
            result['reference_comparison']=dict(source='thesis Table3.1 HS99UU',ADM_energy=.980124,absolute_tolerance=2e-4,energy_error=error,energy_agreement=bool(error<2e-4),exact_historical_reproduction=False,departure='modern superposed-metric trace projection; horizon mass/spin unmeasured')
        elif args.stage=='highboost':
            result['reference_comparison']=dict(source='specified local Gamma=sqrt5 benchmark',exact_historical_reproduction=False,departure='thesis Table4.3 lacks complete bare inputs and uses historical step stuffing')
        if args.stage in ('highspin','highboost'):
            result['accepted_high_regime']=bool(result['passed_strict'] and (args.stage!='highspin' or result['reference_comparison']['energy_agreement']))
    result['stage']=args.stage
    report=json.loads(REPORT.read_text()) if REPORT.exists() else report
    report[label]=result;report['metadata']=dict(library=str(backend.path),library_sha256=backend.library_sha256(),
        numpy_version=np.__version__,python_version=sys.version,cpu_threads=1,date=datetime.now(timezone.utc).isoformat(),horizon_enclosure_verified=False)
    REPORT.write_text(json.dumps(report,indent=2)+'\n')
    accepted=result.get('accepted_high_regime',False) if args.stage in ('highspin','highboost') else result['passed']
    return 0 if accepted else 1

if __name__=='__main__':raise SystemExit(main())
