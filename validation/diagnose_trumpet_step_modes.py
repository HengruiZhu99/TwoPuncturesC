"""Localize a retained failed Newton update; mode removal is diagnostic only."""
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Backend, Config
from checkpoint_export import read_checkpoint
from prolong import for_backend
p=argparse.ArgumentParser(description=__doc__)
for name in ('library','checkpoint','points','final','output'):p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--shape',type=int,nargs=3,required=True)
a=p.parse_args()
if a.output.exists():raise FileExistsError(a.output)
b=Backend(str(a.library.resolve()));c,u,meta=read_checkpoint(a.checkpoint)
if b.parameterization()!=meta['parameterization']:raise ValueError('basis/map mismatch')
old=list(c.n);c.n[:]=a.shape
with np.load(a.points) as z:x=z['xyz'].copy()
with np.load(a.final) as z:final=z['unknowns'].copy()
base=for_backend(b,u,old,a.shape);delta=(final-base).reshape(a.shape[2],a.shape[1],a.shape[0],4)
half=a.shape[2]//2
result=dict(purpose=__doc__,binary_acceptance=False,checkpoint=meta,library_sha256=b.loaded_sha256,
 final_sha256=hashlib.sha256(a.final.read_bytes()).hexdigest(),
 driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),rows=[])
with b.create_sampler(c) as s:
 s.set_unknowns(base);initial=s.sample_with_derivatives(x)
 for mode in range(half+1):
  slots=[mode] if mode in (0,half) else [mode,half+mode]
  d=np.zeros_like(delta);d[slots]=delta[slots]
  s.set_unknowns(base+d.ravel());v=s.sample_with_derivatives(x)
  result['rows'].append(dict(mode=mode,modal_step_linf=np.max(abs(d[slots]),axis=(0,1,2)).tolist(),
    correction_change_linf=np.max(abs(v['correction']-initial['correction']),axis=0).tolist(),
    Kij_change_linf=float(np.max(abs(v['Kij']-initial['Kij']))),
    dgamma_change_linf=float(np.max(abs(v['dgamma']-initial['dgamma'])))))
a.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result['rows'],indent=2))
