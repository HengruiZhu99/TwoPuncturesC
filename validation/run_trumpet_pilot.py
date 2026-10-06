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
from checkpoint_export import write_checkpoint

def run(library,output,n,nphi):
    output.mkdir(parents=True,exist_ok=False)
    b=Backend(str(library.resolve()));c=trumpet_moderate(b,n,nphi)
    result=dict(config={**as_dict(c),'seed_family':c.seed_family},library_sha256=b.library_sha256(),
        kind='moderate_trumpet_diagnostic_pilot',binary_acceptance=False,completed=False,
        criteria=dict(exterior_HM_rms=1e-6,exterior_HM_max=1e-4),
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__).resolve(),ROOT/'examples/trumpet_configs.py',ROOT/'validation/run_validation.py',ROOT/'validation/physical.py']})
    def save(): (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    save();start=time.monotonic()
    with b.create(c,execution='kokkos',geometry='host') as s:
        result['setup_seconds']=time.monotonic()-start;save()
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
        result['linear_history']=s.linear_history();result['work_statistics']=s.work_statistics()
    result.update(completed=True,total_seconds=time.monotonic()-start,process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    save();print(json.dumps({k:result[k] for k in ('diagnostics','physical','minimum_psi','total_seconds')},indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--library',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--n',type=int,default=24);ap.add_argument('--nphi',type=int,default=8);a=ap.parse_args();run(a.library,a.output,a.n,a.nphi)
