"""Native collocation H/M versus independent Cartesian manufactured Hessians.

Tests the full mapped Cartesian second-derivative path, including Fourier
partners and the Nyquist cosine, rather than only coefficient derivatives.
"""
import argparse,json
from pathlib import Path
import numpy as np
from hispid import Backend,Hole
from cartesian_modes import oracle

p=argparse.ArgumentParser();p.add_argument('--library',required=True)
p.add_argument('--levels',default='16,32,64');p.add_argument('--output',required=True)
a=p.parse_args();b=Backend(a.library);records=[];amp=1e-4
for n in map(int,a.levels.split(',')):
    cfg=b.config();np_=16;half=np_//2;cfg.n[:]=[n,n,np_]
    cfg.conformal_choice=0;cfg.inner_flatten=0;cfg.far_radius=0
    cfg.omega[:]=[0,0];cfg.inner_max[:]=[0,0]
    cfg.hole[0]=Hole(.6,(3,0,0));cfg.hole[1]=Hole(.4,(-3,0,0))
    raw=-np.cos(np.pi*(np.arange(n)+.5)/n);sigma=(1+raw)/2
    maps=b.parameterization_maps();lam=maps['radial_stretch'];kap=maps['angular_stretch']
    t=lam*sigma/(1-(1-lam)*sigma);eta=np.tanh(kap*raw)/np.tanh(kap) if kap else raw
    with b.create(cfg) as s:
        baseline=s.equation_samples();x=baseline['xyz']
        baseline_conformal=baseline['physical_equivalent'].copy()
        baseline_conformal[:,0]*=-baseline['psi']**5/8
        baseline_conformal[:,1:]*=baseline['psi'][:,None]**10
        for m in (0,1,2,3,4,5,6,8):
            exponent=m if m<=4 else 3 if m%2 else 4
            for sine in (False,True):
                if sine and m in (0,half):continue
                mode=m if not sine else half+m;normal=np.sqrt((1 if m in (0,half) else 2)/np_)
                P=amp*(t[None,:]*(1-eta[:,None]**2))**((m-exponent)//2)/normal
                v,dv,ddv,seed,_=oracle(x,cfg.hole,3.,amp,m,sine)
                for component in (0,1,2,3):
                    values=np.zeros((np_,n,n,4));values[mode,:,:,component]=P;s.set_unknowns(values.ravel())
                    out=s.equation_samples();psi=seed+v if component==0 else seed
                    H=-8*np.trace(ddv,axis1=1,axis2=2)/psi**5 if component==0 else -(2*np.sum(dv*dv,axis=1)+(2/3)*dv[:,component-1]**2)/psi**12
                    M=np.zeros((len(x),3))
                    if component:
                        M=ddv[:,component-1,:]/3;M[:,component-1]+=np.trace(ddv,axis1=1,axis2=2);M/=psi[:,None]**10
                    expected=np.c_[H,M];difference=out['physical_equivalent']-expected
                    # Convert the discrepancy back to conformal operator units
                    # so large puncture psi cannot hide mapped Hessian errors.
                    conformal=difference.copy();conformal[:,0]*=-psi**5/8;conformal[:,1:]*=psi[:,None]**10
                    # The seed's numerically cancelled Laplacian is a fixed
                    # additive source, independent of the tested correction.
                    # Subtract its separately measured zero-correction value
                    # for this linear mapped-operator check; preserve the full
                    # physical discrepancy and seed floor separately below.
                    conformal-=baseline_conformal
                    scale=expected.copy();scale[:,0]*=-psi**5/8;scale[:,1:]*=psi[:,None]**10
                    normalized=np.max(abs(conformal),axis=0)/np.maximum(np.max(abs(scale),axis=0),amp)
                    row=dict(resolution=list(cfg.n),m=m,sine=sine,component=component,
                        physical_error_max=np.max(abs(difference),axis=0).tolist(),
                        conformal_error_max=np.max(abs(conformal),axis=0).tolist(),
                        normalized_conformal_error=normalized.tolist(),
                        conformal_factor_error=float(np.max(abs(out['psi']-psi))),
                        seed_physical_error_max=np.max(abs(baseline['physical_equivalent']),axis=0).tolist(),
                        seed_conformal_cancellation_max=np.max(abs(baseline_conformal),axis=0).tolist(),
                        passed=bool(np.max(normalized)<1e-8 and np.max(abs(difference))<1e-8))
                    records.append(row);print(json.dumps(row),flush=True)
report=dict(library_sha256=b.library_sha256(),parameterization=b.parameterization(),
    collocation_maps=b.parameterization_maps(),amplitude=amp,records=records,
    all_levels_passed=all(r['passed'] for r in records),
    passed=all(r['passed'] for r in records if r['resolution'][0]==max(v['resolution'][0] for v in records)))
Path(a.output).write_text(json.dumps(report,indent=2)+'\n');raise SystemExit(0 if report['passed'] else 1)
