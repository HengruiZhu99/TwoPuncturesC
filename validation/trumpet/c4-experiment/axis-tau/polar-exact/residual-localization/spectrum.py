"""Meridional DCT energy of retained weighted residuals; no solver evaluation."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--arrays', type=Path, required=True)
p.add_argument('--manifest', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
if a.output.exists():
    raise FileExistsError(a.output)
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
meta = json.loads(a.manifest.read_text())
assert sha(a.arrays) == meta['residual_arrays_sha256']
na, nb, nk = meta['shape']

def spectrum(values, axis):
    # An orthonormal DCT-II on the Chebyshev sampling angles. The sign change
    # from increasing -cos(theta) nodes affects coefficients, not their energy.
    x = np.moveaxis(values, axis, -1)
    n = x.shape[-1]
    mirrored = np.concatenate((x, x[..., ::-1]), axis=-1)
    transformed = np.fft.rfft(mirrored, axis=-1)[..., :n]
    transformed *= np.exp(-1j*np.pi*np.arange(n)/(2*n))
    coeff = transformed.real / np.sqrt(2*n)
    coeff[..., 0] /= np.sqrt(2)
    energy = np.sum(coeff**2, axis=tuple(range(coeff.ndim-1)))
    original = float(np.sum(values**2))
    parseval = abs(float(energy.sum())/original-1)
    assert parseval < 1e-12
    fraction = energy/energy.sum()
    return dict(parseval_relative_error=parseval,
                upper_half_fraction=float(fraction[n//2:].sum()),
                upper_quarter_fraction=float(fraction[3*n//4:].sum()),
                last_eight_fraction=float(fraction[-8:].sum()),
                median_degree=int(np.searchsorted(np.cumsum(fraction), .5)),
                energy_fraction=fraction.tolist())

out = dict(scope='weighted residual DCT-II energy; not a physical constraint norm',
           checkpoint_sha256=meta['checkpoint_sha256'],
           library_sha256=meta['library_sha256'], arrays_sha256=sha(a.arrays),
           driver_sha256=sha(Path(__file__)), shape=meta['shape'])
with np.load(a.arrays, allow_pickle=False) as arrays:
    for name in ('cuda', 'reference'):
        r = arrays[name].reshape(nk, nb, na, 4)
        assert np.isfinite(r).all()
        out[name] = {label: spectrum(r, axis)
                     for label, axis in (('radial', 2), ('polar', 1))}
a.output.write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps({name:{axis:{k:v for k,v in data.items() if k!='energy_fraction'}
                       for axis,data in out[name].items()}
                  for name in ('cuda','reference')}, indent=2))
