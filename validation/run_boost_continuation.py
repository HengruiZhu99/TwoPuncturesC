"""A recorded velocity continuation, followed by unchanged physical gates.

Each stage changes the seed velocity and the Lorentz-contracted g window.
Only the converged fields become an initial guess for the next stage. This
does not change the final benchmark's free data or acceptance thresholds.
"""
import argparse,json,time
import numpy as np
from hispid import Backend
from configs import highboost,as_dict
from checkpoints import ROOT,library_sha,PARAMETERIZATION,restore
from run_validation import solve_case,REPORT,RAW

p=argparse.ArgumentParser();p.add_argument('--library',required=True)
p.add_argument('--levels',default='48:8,80:8,112:8')
p.add_argument('--velocities',default='0,.3,.5,.65,.75,.82,.86,.8944271909999159')
p.add_argument('--target-speed',type=float,default=2/np.sqrt(5))
p.add_argument('--label',default='highboost_continued_far0')
p.add_argument('--restart',type=int,default=80);p.add_argument('--max-krylov',type=int,default=2000)
p.add_argument('--inner-flatten',type=int,choices=(0,1),default=1)
p.add_argument('--memory-mib',type=int,default=2048)
p.add_argument('--inner-min-factor',type=float,default=.2)
p.add_argument('--inner-max-factor',type=float,default=.4)
p.add_argument('--continuation-from',help='saved continuation JSON; start from its last converged velocity and bound checkpoint')
args=p.parse_args();backend=Backend(args.library);sha=library_sha(backend)
report=json.loads(REPORT.read_text())
def gate(label):
    r=report.get(label,{})
    return r.get('passed') and r.get('records') and all(x['library_sha256']==sha for x in r['records'])
moderate_gate=any(gate(label) for label,result in report.items()
                  if isinstance(result,dict) and (result.get('stage')=='moderate' or label=='moderate_far0_stable'))
if not gate('seeds') or not moderate_gate or not gate('covariance') or not gate(report['covariance'].get('source_case','')):
    raise SystemExit('seed, moderate and covariance gates for this library are required')
levels=[tuple(map(int,x.split(':'))) for x in args.levels.split(',')]
velocities=list(map(float,args.velocities.split(',')))
if not 0<=args.target_speed<1 or not velocities or not np.all(np.diff(velocities)>0) or abs(velocities[-1]-args.target_speed)>1e-14:
    raise ValueError('continuation must increase to --target-speed')
if not 0<=args.inner_min_factor<args.inner_max_factor<1:
    raise ValueError('g window factors must satisfy 0<=min<max<1 (seed-radius screening only)')
def factory(backend,n,nphi,q=None):
    if q is None:q=args.target_speed/np.sqrt(1-args.target_speed**2)
    cfg=highboost(backend,n,nphi,q);cfg.far_radius=0
    cfg.inner_flatten=args.inner_flatten;cfg.memory_limit_mib=args.memory_mib
    rh=.25/np.sqrt(1+q*q)
    cfg.inner_min[:]=[args.inner_min_factor*rh]*2;cfg.inner_max[:]=[args.inner_max_factor*rh]*2
    cfg.krylov_restart=args.restart;cfg.max_krylov=args.max_krylov;cfg.max_newton=24
    return cfg
history=[];unknowns=None;last=None;RAW.mkdir(exist_ok=True,parents=True)
if args.continuation_from:
    source=json.loads((ROOT/args.continuation_from).read_text())
    converged=[s for s in source['stages'] if s['diagnostics']['status']==0]
    if not converged:raise ValueError('continuation has no converged checkpoint')
    initial=converged[-1];_,unknowns=restore(backend,initial)
    if abs(velocities[0]-initial['velocity'])>1e-14:
        raise ValueError('initial velocity must match the last converged checkpoint')
elif velocities[0]!=0:raise ValueError('fresh continuation must start at v=0')
for i,v in enumerate(velocities):
    q=v/np.sqrt(1-v*v)
    cfg=factory(backend,*levels[0],q);start=time.monotonic()
    stage=f'{args.label}_velocity{i}'
    with backend.create(cfg) as s:
        if unknowns is not None:s.set_unknowns(unknowns)
        diagnostics=s.solve();unknowns=s.unknowns()
    last=dict(case=stage,resolution=list(cfg.n),config=as_dict(cfg),library_sha256=sha,
              unknown_parameterization=PARAMETERIZATION,diagnostics=diagnostics,
              velocity=v,seconds=time.monotonic()-start)
    np.savez_compressed(RAW/f'{stage}_{cfg.n[0]}_{cfg.n[2]}.npz',unknowns=unknowns)
    history.append(last)
    evidence=dict(label=args.label,library_sha256=sha,stages=history,complete=False)
    (ROOT/f'validation/{args.label}_continuation.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print('velocity',v,diagnostics,flush=True)
    if diagnostics['status']!=0:raise SystemExit('continuation failed; iterate and evidence preserved')
result=solve_case(backend,factory,levels,args.label,horizon_scaled=True,adaptive_steps=True,initial_record=last)
result.update(stage='highboost',continuation_evidence=f'validation/{args.label}_continuation.json',
              accepted_high_regime=bool(result['passed_strict']),
              reference_comparison=dict(source='specified local speed benchmark',target_speed=args.target_speed,exact_historical_reproduction=False))
report=json.loads(REPORT.read_text());report[args.label]=result;REPORT.write_text(json.dumps(report,indent=2)+'\n')
evidence['complete']=True;(ROOT/f'validation/{args.label}_continuation.json').write_text(json.dumps(evidence,indent=2)+'\n')
raise SystemExit(0 if result['accepted_high_regime'] else 1)
