"""Targeted retained-polynomial mode sensitivity; never a solved data product."""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from hispid import Backend
from checkpoint_export import read_checkpoint
from physical import constraints, norms
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--library',required=True);p.add_argument('--pilot',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
if a.output.exists():raise FileExistsError(a.output)
c,u,meta=read_checkpoint(a.pilot/'diagnostic.checkpoint')
b=Backend(a.library)
if b.parameterization()!=meta['parameterization']:raise ValueError('basis mismatch')
if list(c.n)!=[256,512,16]:raise ValueError('this diagnosis is scoped to the retained moderate finest grid')
old=json.loads((a.pilot/'result.json').read_text());start=old['physical']['near']['count'];stop=start+old['physical']['bulk']['count']
with np.load(a.pilot/'physical.npz') as z:x=z['xyz'][start:stop].copy();step=z['step'][start:stop].copy();old_H=z['H'][start:stop].copy()
result=dict(purpose=__doc__,checkpoint=meta,library_sha256=b.library_sha256(),driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),binary_acceptance=False,rows=[])
with b.create_sampler(c) as s:
 for label,removed in [('unchanged',[]),('scalar_m6_removed',[6,14])]:
  values=u.copy().reshape(16,512,256,4)
  for k in removed:values[k,:,:,0]=0
  s.set_unknowns(values.ravel());begin=time.monotonic();r=constraints(s.sample,x,step)
  if not all(np.isfinite(r[k]).all() for k in ('H','M')):raise ValueError('nonfinite observer')
  result['rows'].append(dict(label=label,removed_scalar_rows=removed,bulk=norms(r),H=r['H'].tolist(),Mnorm=r['Mnorm'].tolist(),seconds=time.monotonic()-begin,baseline_H_change_linf=float(np.max(abs(r['H']-old_H)))))
  a.output.write_text(json.dumps(result,indent=2)+'\n');print(label,result['rows'][-1]['bulk'],flush=True)
