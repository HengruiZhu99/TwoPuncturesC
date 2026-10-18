"""Read-only mode localization; filtered fields cannot transfer acceptance."""
import argparse,json,sys,time,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Backend
from checkpoint_export import read_checkpoint
from parallel_observer import ParallelObserver
from physical import constraints,norms
p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
if a.output.exists():raise FileExistsError(a.output)
b=Backend(str(a.library.resolve()));c,u,meta=read_checkpoint(a.run/'diagnostic.checkpoint');na,nb,np_=c.n
with np.load(a.run/'physical.npz') as z:x=z['xyz'][48:66];step=z['step'][48:66];original=z['M'][48:66].copy()
v=u.reshape(np_,nb,na,4);half=np_//2;bands={'all':range(half+1),'none':[],'m0':[0],'m1':[1],'m2to4':[2,3,4],'m5plus':range(5,half+1)}
out=dict(scope=__doc__,binary_acceptance=False,library_sha256=b.library_sha256(),checkpoint=meta,records={})
start=time.monotonic()
with b.create_sampler(c) as s:
 for label,ms in bands.items():
  w=v.copy();w[:,:,:,1:]=0
  for m in ms:
   for k in ([m] if m in (0,half) else [m,half+m]):w[k,:,:,1:]=v[k,:,:,1:]
  s.set_unknowns(w.ravel())
  with ParallelObserver(s,16) as sample:q=constraints(sample,x,step)
  out['records'][label]=dict(norms=norms(q),M=q['M'].tolist())
  if label=='all':out['baseline_component_difference_linf']=float(np.max(abs(q['M']-original)))
  out['seconds']=time.monotonic()-start;a.output.write_text(json.dumps(out,indent=2)+'\n')
print({k:r['norms']['M_rms'] for k,r in out['records'].items()})
