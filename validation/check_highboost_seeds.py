"""Independent exact-seed controls at the local Gamma=sqrt5 benchmark boost."""
import argparse,json,time
import numpy as np
from hispid import Backend
from configs import highboost,as_dict
from run_validation import points
from physical import constraints,norms,charges,extrapolate
from checkpoints import ROOT,library_sha

p=argparse.ArgumentParser();p.add_argument('--library',required=True)
p.add_argument('--output',default='validation/highboost_seed_controls.json')
args=p.parse_args();backend=Backend(args.library);cfg=highboost(backend)
xyz,near,bulk=points(cfg,True)
output=dict(library_sha256=library_sha(backend),config=as_dict(cfg),records=[],passed=False)
for i,hole in enumerate(cfg.hole):
    x=np.r_[xyz[i*(near//2):(i+1)*(near//2)],xyz[near:near+bulk]]
    step=np.minimum(.002,.001*np.linalg.norm(x-np.array(hole.center),axis=1))
    sample=lambda x,h=hole:backend.seed(h,x,cfg.conformal_choice)
    sequence=[];start=time.monotonic()
    for factor in (2,1,.5):
        r=constraints(sample,x,factor*step)
        sequence.append(dict(step_factor=factor,norms=norms(r),H=r['H'].tolist(),M=r['M'].tolist()))
    radii=[40.,80.,160.,320.,640.,1280.]
    v=np.array(hole.velocity);G=1/np.sqrt(1-v@v)
    expected=np.r_[G*hole.mass,G*hole.mass*v,[0,0,0]]
    quadratures=[]
    for nt,np_ in ((16,32),(24,48),(32,64)):
        q=[charges(sample,R,center=list(hole.center),ntheta=nt,nphi=np_) for R in radii]
        fit=extrapolate(radii[-4:],q[-4:]);error=abs(fit-expected)
        quadratures.append(dict(ntheta=nt,nphi=np_,charges=[v.tolist() for v in q],
            extrapolated=fit.tolist(),charge_error=error.tolist()))
    change=max(abs(np.array(quadratures[-1]['extrapolated'])-quadratures[-2]['extrapolated']))
    radial=[dict(fit_radii=radii[i:i+4],extrapolated=extrapolate(radii[i:i+4],q[i:i+4]).tolist()) for i in range(3)]
    radial_change=max(abs(np.array(radial[-1]['extrapolated'])-radial[-2]['extrapolated']))
    passed=sequence[-1]['norms']['H_rms']<1e-7 and sequence[-1]['norms']['M_rms']<1e-7 and max(error)<2e-5 and change<2e-5 and radial_change<2e-5
    output['records'].append(dict(hole=as_dict(hole),points=x.tolist(),verifier_steps=step.tolist(),sequence=sequence,
        radii=radii,charges=[v.tolist() for v in q],extrapolated=fit.tolist(),expected=expected.tolist(),
        charge_error=error.tolist(),fit_radii=radii[-4:],radial_sequence=radial,last_radial_change=float(radial_change),
        quadrature_sequence=quadratures,last_quadrature_change=float(change),
        seconds=time.monotonic()-start,passed=bool(passed)))
    print('boosted seed',i,sequence[-1]['norms'],'charge error',error,flush=True)
    (ROOT/args.output).write_text(json.dumps(output,indent=2)+'\n')
output['passed']=all(r['passed'] for r in output['records'])
(ROOT/args.output).write_text(json.dumps(output,indent=2)+'\n')
raise SystemExit(0 if output['passed'] else 1)
