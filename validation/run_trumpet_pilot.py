"""One diagnosed moderate binary pilot on the selected production path.

Always exports diagnostic data; a single grid cannot establish convergence.
"""
import argparse,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'examples')]
from hispid import Backend
from configs import as_dict
from trumpet_configs import trumpet_moderate
from run_validation import points
from physical import constraints,norms
from checkpoint_export import write_checkpoint,read_checkpoint
from prolong import for_backend

def run(library,output,n,nphi,initial=None,memory_mib=2048,npolar=None):
    output.mkdir(parents=True,exist_ok=False)
    b=Backend(str(library.resolve()));c=trumpet_moderate(b,n,nphi);c.memory_limit_mib=memory_mib
    if npolar is not None:c.n[1]=npolar
    if any(k<4 or k>256 for k in c.n) or c.n[2]%2:raise ValueError("native grid dimensions must be4..256 with even nphi")
    result=dict(config={**as_dict(c),'seed_family':c.seed_family},library_sha256=b.library_sha256(),
        kind='moderate_trumpet_diagnostic_pilot',binary_acceptance=False,completed=False,
        criteria=dict(exterior_HM_rms=1e-6,exterior_HM_max=1e-4),
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__).resolve(),ROOT/'examples/trumpet_configs.py',ROOT/'validation/run_validation.py',ROOT/'validation/physical.py']})
    def save(): (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    save();start=time.monotonic()
    with b.create(c,execution='kokkos',geometry='host') as s:
        result['setup_seconds']=time.monotonic()-start;save()
        if initial:
            old,values,meta=read_checkpoint(initial)
            if meta['source_library_sha256']!=b.loaded_sha256 or old.seed_family!=c.seed_family or meta['parameterization']!=b.parameterization():raise ValueError('incompatible initial checkpoint')
            excluded={'n','tolerance','max_newton','max_krylov','krylov_restart','memory_limit_mib'}
            if any(as_dict(old)[k]!=v for k,v in as_dict(c).items() if k not in excluded):raise ValueError('initial checkpoint has different physical free data')
            s.set_unknowns(for_backend(b,values,list(old.n),list(c.n)))
            result['initial_checkpoint']={**meta,'acceptance_inherited':False};save()
        result['diagnostics']=s.solve(krylov='gmres',linear_rtol=.1)
        values=s.unknowns();np.savez_compressed(output/'solve.npz',unknowns=values)
        result['solve_artifact_sha256']=hashlib.sha256((output/'solve.npz').read_bytes()).hexdigest();save()
        result['checkpoint']=write_checkpoint(output/'diagnostic.checkpoint',c,values,b.library_sha256(),'diagnostic',b.parameterization())
        x,near,bulk=points(c,True);regions=dict(near=np.arange(len(x))<near,bulk=(np.arange(len(x))>=near)&(np.arange(len(x))<near+bulk),modified=np.arange(len(x))>=near+bulk)
        distance=np.min([np.linalg.norm(x-np.array(h.center),axis=1) for h in c.hole],axis=0)
        step=np.minimum(.002,.002*distance)
        physical=constraints(s.sample,x,step);sampled=s.sample(x)
        result['physical']={key:norms(physical,mask) for key,mask in regions.items()}
        result['exterior_stencil_unmodified']=bool(np.all(physical['stencil_attenuation_all_one'][:near+bulk]))
        result['minimum_psi']=float(np.min(sampled['psi']));result['minimum_metric_eigenvalue']=float(np.min(physical['min_metric_eigenvalue']))
        np.savez_compressed(output/'physical.npz',xyz=x,step=step,**physical)
        radii=[256.,512.,1024.]
        result['charge_radii']=radii
        result['charges']=[s.charges(r,ntheta=2*c.n[1],nphi=32).tolist() for r in radii]
        from physical import extrapolate
        result['charges_extrapolated']=extrapolate(radii,result['charges']).tolist()
        result['linear_history']=s.linear_history();result['work_statistics']=s.work_statistics()
    result.update(completed=True,total_seconds=time.monotonic()-start,process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    save();print(json.dumps({k:result[k] for k in ('diagnostics','physical','minimum_psi','total_seconds')},indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--library',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--n',type=int,default=24);ap.add_argument('--nphi',type=int,default=8);ap.add_argument('--initial',type=Path);ap.add_argument('--memory-mib',type=int,default=2048);ap.add_argument('--npolar',type=int);a=ap.parse_args();run(a.library,a.output,a.n,a.nphi,a.initial,a.memory_mib,a.npolar)
