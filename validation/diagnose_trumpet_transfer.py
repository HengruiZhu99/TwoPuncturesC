"""Compare retained fields before/after modal prolongation; no elliptic solve."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'python'), str(ROOT/'validation')]
from hispid import Backend, Config
from checkpoint_export import read_checkpoint
from prolong import for_backend

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--library', type=Path, required=True)
p.add_argument('--checkpoint', type=Path, required=True)
p.add_argument('--points', type=Path, required=True)
p.add_argument('--shape', type=int, nargs=3, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
if a.output.exists(): raise FileExistsError(a.output)
b = Backend(str(a.library.resolve()))
c, u, meta = read_checkpoint(a.checkpoint)
if meta['parameterization'] != b.parameterization(): raise ValueError('basis/map mismatch')
with np.load(a.points) as z: xyz = z['xyz'].copy()
fine = Config.from_buffer_copy(c)
fine.seed_family = c.seed_family
fine.n[:] = a.shape
v = for_backend(b, u, list(c.n), a.shape)
fields = []
start = time.monotonic()
for config, values in ((c,u),(fine,v)):
    with b.create_sampler(config) as s:
        s.set_unknowns(values)
        fields.append(s.sample_with_derivatives(xyz))
errors = {}
for k, x in fields[0].items():
    y = fields[1][k]
    if not np.isfinite(x).all() or not np.isfinite(y).all(): raise ValueError('nonfinite '+k)
    errors[k] = dict(absolute_linf=float(np.max(abs(x-y))),
                    scaled_linf=float(np.max(abs(x-y)/(1+abs(x)))))
r = dict(purpose=__doc__, checkpoint=meta, library_sha256=b.loaded_sha256,
         driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         points_sha256=hashlib.sha256(a.points.read_bytes()).hexdigest(),
         old_shape=list(c.n), new_shape=a.shape, point_count=len(xyz),
         binary_acceptance=False, errors=errors, seconds=time.monotonic()-start)
a.output.write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
