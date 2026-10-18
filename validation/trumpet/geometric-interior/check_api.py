"""Separate-process QI preservation witness and one changed-interior JVP control."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'examples')]
from hispid import Backend
from trumpet_configs import trumpet_moderate_wide
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--baseline',type=Path);a=p.parse_args()
if a.output.exists():raise FileExistsError(a.output)
b=Backend(str(a.library.resolve()));c=trumpet_moderate_wide(b,12,16);c.seed_family='qi'
with b.create(c) as s:
 t=np.arange(s.size);u=1e-4*np.sin(.13*t);v=1e-3*np.cos(.27*t)
 r,j=s.residual(u),s.jvp(u,v)
np.savez_compressed(a.output.with_suffix('.npz'),residual=r,jvp=j)
q={'library_sha256':b.loaded_sha256}
if a.baseline:
 with np.load(a.baseline.with_suffix('.npz')) as old:
  q.update(qi_residual_bitwise_equal=bool(np.array_equal(old['residual'],r)),qi_jvp_bitwise_equal=bool(np.array_equal(old['jvp'],j)))
 assert q['qi_residual_bitwise_equal'] and q['qi_jvp_bitwise_equal']
 c.seed_family='trumpet_r0_m'
 with b.create(c) as s:
  j=s.jvp(u,v)
  q['centered_jvp_relative_errors']=[float(np.linalg.norm((s.residual(u+h*v)-s.residual(u-h*v))/(2*h)-j)/np.linalg.norm(j)) for h in [.01,.005]]
  assert max(q['centered_jvp_relative_errors'])<1e-7
q['passed']=True
a.output.write_text(json.dumps(q,indent=2)+'\n');print(q)
