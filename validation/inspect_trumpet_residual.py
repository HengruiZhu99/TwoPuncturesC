"""Locate a retained checkpoint's weighted residual without taking a solve step."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from hispid import Backend
from checkpoint_export import read_checkpoint

p = argparse.ArgumentParser(description=__doc__)
for name in ('library', 'source-library', 'checkpoint', 'output'):
    p.add_argument('--' + name, type=Path, required=True)
p.add_argument('--compare-reference', action='store_true')
a = p.parse_args()
partial_json = a.output.with_suffix('.cuda.json')
partial_array = a.output.with_suffix('.cuda.npy')
if any(path.exists() for path in (a.output, a.output.with_suffix('.npz'),
                                  partial_json, partial_array)):
    raise FileExistsError(a.output)
c, values, meta = read_checkpoint(a.checkpoint)
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
driver_sha = sha(Path(__file__))
started = time.monotonic()
def progress(stage):
    print(json.dumps(dict(stage=stage, elapsed_seconds=time.monotonic()-started)),
          flush=True)
assert sha(a.source_library) == meta['source_library_sha256']
b = Backend(str(a.library.resolve()))
assert b.parameterization() == meta['parameterization']
progress('cuda_context_start')
with b.create(c, execution='kokkos', geometry='host') as s:
    progress('cuda_context_ready')
    s.set_unknowns(values)
    residual = s.residual(values)
progress('cuda_residual_complete')
assert np.isfinite(residual).all()
na, nb, np_ = map(int, c.n)
r = residual.reshape(np_, nb, na, 4)
energy = np.sum(r*r, axis=3)
total = float(energy.sum())
assert total > 0
top = np.argsort(energy.ravel())[-20:][::-1]
rows = []
for q in top:
    k, j, i = map(int, np.unravel_index(q, energy.shape))
    rows.append(dict(i=i, j=j, k=k, residual=r[k,j,i].tolist(),
                     fraction=float(energy[k,j,i]/total)))
layers = {}
for width in (1, 2, 4, 8):
    radial = np.minimum(np.arange(na), na-1-np.arange(na)) < width
    polar = np.minimum(np.arange(nb), nb-1-np.arange(nb)) < width
    layers[str(width)] = dict(radial=float(energy[:,:,radial].sum()/total),
                             polar=float(energy[:,polar,:].sum()/total))
out = dict(library_sha256=b.loaded_sha256, checkpoint_sha256=sha(a.checkpoint),
           source_library_sha256=meta['source_library_sha256'], shape=list(c.n),
           parameterization=b.parameterization(), residual_l2=np.sqrt(total),
           residual_linf=float(abs(r).max()), largest_points=rows,
           component_energy_fraction=(np.sum(r*r,axis=(0,1,2))/total).tolist(),
           radial_energy_fraction=(energy.sum(axis=(0,1))/total).tolist(),
           polar_energy_fraction=(energy.sum(axis=(0,2))/total).tolist(),
           azimuthal_point_energy_fraction=(energy.sum(axis=(1,2))/total).tolist(),
           endpoint_layers=layers, driver_sha256=driver_sha,
           scope='weighted collocation residual localization only', binary_acceptance=False)
if a.compare_reference:
    # Keep a completed stage if the slower serial reference exceeds wall time.
    # These are explicitly partial artifacts, never a successful comparison.
    np.save(partial_array, residual, allow_pickle=False)
    partial = dict(out, stage='cuda_only_reference_pending',
                   residual_array_sha256=sha(partial_array))
    partial_json.write_text(json.dumps(partial, indent=2) + '\n')
    # Sequential contexts avoid retaining both full geometry caches at once.
    # Neither execution path is treated as an exact continuum reference.
    progress('reference_context_start')
    with b.create(c, execution='reference', geometry='host') as s:
        progress('reference_context_ready')
        s.set_unknowns(values)
        reference = s.residual(values)
    progress('reference_residual_complete')
    assert np.isfinite(reference).all()
    delta = reference-residual
    out['reference_comparison'] = dict(
        reference_l2=float(np.linalg.norm(reference)),
        reference_linf=float(np.max(abs(reference))),
        difference_l2=float(np.linalg.norm(delta)),
        difference_linf=float(np.max(abs(delta))),
        difference_relative_to_cuda=float(np.linalg.norm(delta)/np.sqrt(total)),
        component_difference_l2=np.linalg.norm(delta.reshape(-1,4),axis=0).tolist(),
        unknown_component_linf=np.max(abs(values.reshape(-1,4)),axis=0).tolist())
    np.savez_compressed(a.output.with_suffix('.npz'), cuda=residual, reference=reference)
    out['residual_arrays_sha256'] = sha(a.output.with_suffix('.npz'))
a.output.write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({k:out[k] for k in ('shape','residual_l2','residual_linf','component_energy_fraction','endpoint_layers')},indent=2))
