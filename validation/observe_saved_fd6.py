"""Higher-accuracy physical observers of immutable saved fields; no solve."""
import argparse,ctypes as C,copy,hashlib,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from physical import constraints,norms
from parallel_observer import ParallelObserver
p=argparse.ArgumentParser();p.add_argument('--system',choices=('hispid','by'),required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=False);r=json.loads((a.run/'result.json').read_text());assert r.get('completed',True)
with np.load(a.run/'physical.npz') as z:x=z['xyz'][:66].copy();step=z['step'][:66].copy();original={k:z[k][:66].copy() for k in ('H','M','Mnorm','R','trace','K2')}
assert hashlib.sha256((a.run/'physical.npz').read_bytes()).hexdigest()==r.get('observer_artifact_sha256',r.get('physical_sha256'))
result=dict(system=a.system,scope=__doc__,source_result_sha256=hashlib.sha256((a.run/'result.json').read_bytes()).hexdigest(),method='FD6, four times the stored FD4 step; same66 exterior points',nonlinear_solve_repeated=False,source_stopping_gate=r['diagnostics']['converged'])
start=time.monotonic()
def save():(a.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
save()
def observe(sample):
 old=constraints(sample,x,step);errors={k:float(np.max(abs(old[k]-original[k]))) for k in original};result['fd4_replay_errors']=errors;save()
 # The same native image and coefficients must reproduce the old observer.
 assert max(errors.values())<1e-13,errors
 new=constraints(sample,x,4*step,order=6);np.savez_compressed(a.output/'physical.npz',xyz=x,step=4*step,**new)
 result['physical']={region:norms(new,mask) for region,mask in [('near',np.arange(66)<48),('bulk',np.arange(66)>=48)]}
 result['exterior_unmodified']=bool(np.all(new['stencil_attenuation_all_one']))
if a.system=='hispid':
 from hispid import Backend
 from checkpoint_export import read_checkpoint
 c,u,meta=read_checkpoint(a.run/'diagnostic.checkpoint');b=Backend(str(a.library.resolve()));assert b.library_sha256()==r['library_sha256']
 result['checkpoint']=meta;result['library_sha256']=b.library_sha256();result['resolution']=list(c.n)
 with b.create_sampler(c) as s:
  s.set_unknowns(u)
  with ParallelObserver(s,16) as sample:observe(sample)
else:
 from bowen_york import Backend,InitialData,PTR
 b=Backend(str(a.library.resolve()));assert b.library_sha256()==r['library_sha256'];cfg=copy.deepcopy(r['config']);cfg['integer']['Newton_maxit']=0;cfg['integer']['verbose']=0
 result['library_sha256']=b.library_sha256();result['resolution']=r['resolution'];result['snapshot_sha256']=hashlib.sha256((a.run/'unknowns.npy').read_bytes()).hexdigest()
 with b.create(cfg) as s:
  s.solve(linear_rtol=.001,krylov='gmres',preconditioner='modal',max_krylov=512,restart=64,execution='reference')
  result['reconstruction_newton_steps']=s.work_statistics()['newton_iterations'];assert result['reconstruction_newton_steps']==0
  v=np.ascontiguousarray(np.load(a.run/'unknowns.npy'),dtype=float);data=C.cast(s.context,C.POINTER(InitialData)).contents;assert v.size==data.ntotal
  C.memmove(data.v.contents.d0,v.ctypes.data,v.nbytes)
  b.lib.SpecCoef.argtypes=[C.c_int]*4+[PTR,PTR];b.lib.SpecCoef.restype=None
  b.lib.SpecCoef(*r['resolution'],1,data.v.contents.d0,data.cf_v.contents.d0)
  assert np.array_equal(s.unknowns(),v)
  observe(s.sample)
result['physical_sha256']=hashlib.sha256((a.output/'physical.npz').read_bytes()).hexdigest();result['seconds']=time.monotonic()-start
(a.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print({k:q['H_rms'] for k,q in result['physical'].items()})
