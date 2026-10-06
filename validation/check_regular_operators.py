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
p.add_argument('--npolar',type=int);p.add_argument('--nphi',type=int,default=16)
p.add_argument('--modes',default='0,1,2,3,4,5,6,8');p.add_argument('--components',default='0,1,2,3')
p.add_argument('--memory-mib',type=int)
p.add_argument('--regularity-cap',type=int,choices=[4,6],default=4)
p.add_argument('--execution',choices=['reference','kokkos'],default='reference')
a=p.parse_args()
if Path(a.output).exists():raise FileExistsError(a.output)
b=Backend(a.library);records=[];amp=1e-4
if not b.parameterization().startswith('modal_P_C'+str(a.regularity_cap-2)+'prolate_'):raise ValueError('declared regularity cap does not match library basis')
for n in map(int,a.levels.split(',')):
    cfg=b.config();np_=a.nphi;nb=a.npolar or n;half=np_//2;cfg.n[:]=[n,nb,np_]
    if a.memory_mib is not None:cfg.memory_limit_mib=a.memory_mib
    modes=list(map(int,a.modes.split(',')));components=list(map(int,a.components.split(',')))
    if any(m<0 or m>half for m in modes) or any(c not in (0,1,2,3) for c in components):raise ValueError('invalid mode/component')
    cfg.conformal_choice=0;cfg.inner_flatten=0;cfg.far_radius=0
    cfg.omega[:]=[0,0];cfg.inner_max[:]=[0,0]
    cfg.hole[0]=Hole(.6,(3,0,0));cfg.hole[1]=Hole(.4,(-3,0,0))
    raw=-np.cos(np.pi*(np.arange(n)+.5)/n);sigma=(1+raw)/2
    maps=b.parameterization_maps();lam=maps['radial_stretch'];kap=maps['angular_stretch']
    t=lam*sigma/(1-(1-lam)*sigma);polar_raw=-np.cos(np.pi*(np.arange(nb)+.5)/nb)
    eta=np.tanh(kap*polar_raw)/np.tanh(kap) if kap else polar_raw
    with b.create(cfg,execution=a.execution) as s:
        baseline=s.equation_samples();x=baseline['xyz']
        baseline_conformal=baseline['physical_equivalent'].copy()
        baseline_conformal[:,0]*=-baseline['psi']**5/8
        baseline_conformal[:,1:]*=baseline['psi'][:,None]**10
        for m in modes:
            exponent=m if m<=a.regularity_cap else a.regularity_cap-1 if m%2 else a.regularity_cap
            for sine in (False,True):
                if sine and m in (0,half):continue
                mode=m if not sine else half+m;normal=np.sqrt((1 if m in (0,half) else 2)/np_)
                P=amp*(t[None,:]*(1-eta[:,None]**2))**((m-exponent)//2)/normal
                v,dv,ddv,seed,_=oracle(x,cfg.hole,3.,amp,m,sine)
                for component in components:
                    values=np.zeros((np_,nb,n,4));values[mode,:,:,component]=P;s.set_unknowns(values.ravel())
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
report=dict(library_sha256=b.library_sha256(),parameterization=b.parameterization(),execution=a.execution,regularity_cap=a.regularity_cap,
    collocation_maps=b.parameterization_maps(),amplitude=amp,records=records,
    all_levels_passed=all(r['passed'] for r in records),
    passed=all(r['passed'] for r in records if r['resolution'][0]==max(v['resolution'][0] for v in records)))
if a.regularity_cap==6:
    import tempfile
    from checkpoint_export import write_checkpoint,read_checkpoint
    from prolong import for_backend
    from remapped_guess import regularity_cap
    small=b.config();small.n[:]=[4,4,16];small.seed_family='trumpet_r0_m'
    source=np.zeros((16,4,4,4));source[6,:,:,0]=1
    with tempfile.TemporaryDirectory() as tmp:
        path=Path(tmp)/'c4.checkpoint'
        write_checkpoint(path,small,source.ravel(),b.library_sha256(),'diagnostic',b.parameterization())
        _,restored,metadata=read_checkpoint(path)
    target=for_backend(b,source.ravel(),[4,4,16],[6,6,16]).reshape(16,6,6,4)
    expected=np.zeros_like(target);expected[6,:,:,0]=1
    error=float(np.max(abs(target-expected)))
    report['metadata']=dict(roundtrip_exact=bool(np.array_equal(restored,source.ravel())),
        parameterization_preserved=metadata['parameterization']==b.parameterization(),
        distinct_from_default=regularity_cap(b.parameterization())!=regularity_cap('modal_P_C2prolate_mapped_v2'),
        constant_mode_prolongation_error=error)
    report['passed']=bool(report['passed'] and report['metadata']['roundtrip_exact'] and
        report['metadata']['parameterization_preserved'] and report['metadata']['distinct_from_default'] and error<1e-14)
Path(a.output).write_text(json.dumps(report,indent=2)+'\n');raise SystemExit(0 if report['passed'] else 1)
