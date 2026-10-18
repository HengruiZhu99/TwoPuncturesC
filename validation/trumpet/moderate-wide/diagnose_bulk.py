"""Saved-field bulk sensitivity only; filtered fields are never solved/accepted data."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Backend
from checkpoint_export import read_checkpoint
from physical import constraints,norms
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--pilot',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--source-library',type=Path);p.add_argument('--full-only',action='store_true');a=p.parse_args()
if a.output.exists():raise FileExistsError(a.output)
b=Backend(str(a.library.resolve()));c,u,meta=read_checkpoint(a.pilot/'diagnostic.checkpoint')
assert (hashlib.sha256(a.source_library.read_bytes()).hexdigest() if a.source_library else b.loaded_sha256)==meta['source_library_sha256']
assert b.parameterization()==meta['parameterization']
old=json.loads((a.pilot/'result.json').read_text());n=old['physical']['near']['count'];nb=old['physical']['bulk']['count']
with np.load(a.pilot/'physical.npz') as raw:x=raw['xyz'][n:n+nb].copy();step=raw['step'][n:n+nb].copy()
na,npol,nphi=map(int,c.n);half=nphi//2;original=u.reshape(nphi,npol,na,4)
rows=[];start=time.monotonic()
with b.create_sampler(c) as s:
 for label,removed,h in ([('full',[],1)] if a.full_only else [('full',[],1),('full_half_step',[],.5)]+[(f'remove_m{m}',[m]+([m+half] if m!=half else []),1) for m in (5,6,7,8)]):
  v=original.copy();v[removed]=0;s.set_unknowns(v.ravel());r=constraints(s.sample,x,step*h)
  row=dict(label=label,removed_modes=removed,bulk=norms(r),H=r['H'].tolist(),Mnorm=r['Mnorm'].tolist())
  rows.append(row);print(label,row['bulk'],flush=True)
  a.output.write_text(json.dumps(dict(scope=__doc__,checkpoint=meta,rows=rows,seconds=time.monotonic()-start,binary_acceptance=False,driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
