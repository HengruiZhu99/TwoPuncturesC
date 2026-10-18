"""Check that the candidate changes interior equations only; no acceptance claim."""
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'examples')]
from hispid import Backend
from trumpet_configs import trumpet_moderate, trumpet_moderate_wide
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
if a.output.exists():raise FileExistsError(a.output)
b=Backend(str(a.library.resolve()));rows=[]
for factory in (trumpet_moderate,trumpet_moderate_wide):
 c=factory(b,24,16)
 with b.create(c) as s:
  idx=np.arange(s.size);u=1e-8*np.sin(.13*idx);v=1e-7*np.cos(.27*idx)
  s.set_unknowns(u);q=s.equation_samples()
  rows.append((q,s.residual(u).reshape(-1,4),s.jvp(u,v).reshape(-1,4)))
x=rows[0][0]['xyz'];centers=np.array([h.center for h in c.hole])
dist=np.linalg.norm(x[:,None,:]-centers[None,:,:],axis=2)
exterior=np.all(dist>=np.array(c.inner_max)[None,:],axis=1)
changed=(rows[0][0]['attenuation']!=rows[1][0]['attenuation'])
assert exterior.sum()>100 and changed.any()
assert np.all(rows[0][0]['attenuation'][exterior]==1) and np.all(rows[1][0]['attenuation'][exterior]==1)
out={'scope':'frozen-field exterior equation identity; no solve or binary acceptance',
     'library_sha256':b.loaded_sha256,'parameterization':b.parameterization(),
     'exterior_points':int(exterior.sum()),'changed_interior_points':int(changed.sum()),
     'inner_min':list(c.inner_min),'inner_max':list(c.inner_max),
     'exterior_residual_bitwise_equal':bool(np.array_equal(rows[0][1][exterior],rows[1][1][exterior])),
     'exterior_jvp_bitwise_equal':bool(np.array_equal(rows[0][2][exterior],rows[1][2][exterior])),
     'interior_residual_changes':bool(np.any(rows[0][1][changed]!=rows[1][1][changed])),
     'driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
assert out['exterior_residual_bitwise_equal'] and out['exterior_jvp_bitwise_equal'] and out['interior_residual_changes']
a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
