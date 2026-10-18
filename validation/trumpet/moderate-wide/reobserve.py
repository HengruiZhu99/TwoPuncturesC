"""Reobserve the unchanged pilot on the declared exterior shells; retain old data."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation'),str(ROOT/'examples')]
from hispid import Backend
from checkpoint_export import read_checkpoint
from run_validation import points
from physical import constraints,norms
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--pilot',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
a.output.mkdir(exist_ok=False,parents=True)
r=json.loads((a.pilot/'result.json').read_text());c,u,meta=read_checkpoint(a.pilot/'diagnostic.checkpoint');b=Backend(str(a.library.resolve()))
assert meta['source_library_sha256']==b.loaded_sha256 and r['diagnostics']['converged'] and r['completed']
# Retain all original scientific settings and solve statistics. This is a new
# observer of the same immutable checkpoint, not a new nonlinear solution.
x,near,bulk=points(c,True,transition_aware=True)
distance=np.min(np.linalg.norm(x[:,None,:]-np.array([h.center for h in c.hole])[None,:,:],axis=2),axis=1);step=np.minimum(.002,.002*distance)
start=time.monotonic()
with b.create_sampler(c) as s:
 s.set_unknowns(u);v=constraints(s.sample,x,step);sample=s.sample(x)
r['physical']={k:norms(v,m) for k,m in dict(near=np.arange(len(x))<near,bulk=(np.arange(len(x))>=near)&(np.arange(len(x))<near+bulk),modified=np.arange(len(x))>=near+bulk).items()}
r['exterior_stencil_unmodified']=bool(np.all(v['stencil_attenuation_all_one'][:near+bulk]))
assert r['exterior_stencil_unmodified']
r['minimum_psi']=float(np.min(sample['psi']));r['minimum_metric_eigenvalue']=float(np.min(v['min_metric_eigenvalue']))
r['reobservation']={'original_result':str(a.pilot/'result.json'),'original_sha256':hashlib.sha256((a.pilot/'result.json').read_bytes()).hexdigest(),'checkpoint_sha256':meta['file_sha256'],'profile':'transition-aware exterior shells: max(legacy shell,1.1 inner_max)','seconds':time.monotonic()-start,'driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
np.savez_compressed(a.output/'physical.npz',xyz=x,step=step,**v)
# Keep a convenience link for continuation; the bound original is not rewritten.
(a.output/'diagnostic.checkpoint').symlink_to((a.pilot/'diagnostic.checkpoint').resolve())
(a.output/'result.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r['physical'],indent=2),flush=True)
