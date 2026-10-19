"""Replay saved binary unknowns in the original paper's solution norm.

No nonlinear solve: compare u/(A-1) on a fixed 6^3 compactified grid,
including infinity, both foci and the axis. Report absolute differences for
the vector components and wherever a relative denominator could be zero.
"""
import argparse
import ctypes as C
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'python'), str(ROOT/'validation')]
from bowen_york import Backend, PTR, ptr
from hispid import Backend as HiSpID
from prolong import chebyshev_evaluation_matrix


def modal_U(values, shape, coords):
    a, B, phi = coords.T
    nr, nb, np_ = shape
    t = a*a
    eta = -2*B/(1+B*B)
    sigma = t/(.2+.8*t)
    zeta = np.arctanh(np.tanh(2.)*eta)/2
    R = chebyshev_evaluation_matrix(nr, 2*sigma-1)
    T = chebyshev_evaluation_matrix(nb, zeta)
    P = np.einsum('qi,qj,kjiv->qkv', R, T, values.reshape(np_, nb, nr, 4), optimize=True)
    answer = np.zeros((len(a), 4))
    for k in range(np_):
        m = k if k <= np_//2 else k-np_//2
        exponent = m if m <= 4 else (3 if m % 2 else 4)
        normal = np.sqrt((1. if m in (0, np_//2) else 2.)/np_)
        basis = normal*(np.cos(m*phi) if k <= np_//2 else np.sin(m*phi))
        factor = 2*(1+a)*(a*np.sqrt(np.maximum(0, 1-eta*eta)))**exponent
        answer += (factor*basis)[:, None]*P[:, k]
    return answer


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study', type=Path, required=True)
    p.add_argument('--by-library', type=Path, required=True)
    p.add_argument('--hispid-library', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    report = json.loads((a.study/'result.json').read_text())
    # Load HiSpID and verify its dependency before opening the explicit BY
    # image. Both archived builds share a dependency SONAME, so loading BY
    # first can substitute its image when the dynamic linker opens HiSpID.
    hb = HiSpID(str(a.hispid_library.resolve()))
    backend = Backend(str(a.by_library.resolve()))
    if backend.library_sha256() != report['libraries']['bowen_york']:
        raise ValueError('source native image differs')
    if hb.library_sha256() != report['libraries']['hispid'] or hb.parameterization_maps() != dict(radial_stretch=.2, angular_stretch=2.) or not hb.parameterization().startswith('modal_P_C2tauC4_focus_map_v9_'):
        raise ValueError('this replay requires the bound C1 .2/polar2 image')
    lib = backend.lib
    lib.SpecCoef.argtypes = [C.c_int]*4+[PTR, PTR]
    lib.SpecCoef.restype = None
    lib.PunctEvalAtArbitPositionFast.argtypes = [PTR, C.c_int]+[C.c_double]*3+[C.c_int]*4
    lib.PunctEvalAtArbitPositionFast.restype = C.c_double
    coords = np.array(np.meshgrid(np.linspace(0, 1, 6), np.linspace(-1, 1, 6), np.arange(6)*2*np.pi/6, indexing='ij')).reshape(3, -1).T
    stored = {}
    hashes = {}
    for rec in report['records']:
        system, shape = rec['system'], rec['n']
        path = a.study/f'{system}{shape[0]}'/'unknowns.npy'
        values = np.ascontiguousarray(np.load(path))
        hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        if system == 'hispid':
            # This replay explicitly supports only the study's unchanged
            # original Mobius .2/polar2, C2-modal basis with W=0.
            if rec['config']['far_radius'] != 0:
                raise ValueError('W!=0 is outside the replay definition')
            U = modal_U(values, shape, coords)
        else:
            coefficients = np.empty_like(values)
            lib.SpecCoef(*shape, 1, ptr(values), ptr(coefficients))
            U = np.array([2*lib.PunctEvalAtArbitPositionFast(ptr(coefficients), 0, 2*x-1, y, phi, 1, *shape) for x, y, phi in coords])[:, None]
        stored[system, shape[0]] = U
    result = dict(method=__doc__, study_sha256=hashlib.sha256((a.study/'result.json').read_bytes()).hexdigest(),
                  library_sha256=backend.library_sha256(), unknown_sha256=hashes,
                  hispid_basis=dict(library_sha256=hb.library_sha256(), parameterization=hb.parameterization(), maps=hb.parameterization_maps()),
                  reference_n=max(report['levels']), records=[])
    for rec in report['records']:
        system, n = rec['system'], rec['n'][0]
        if n == result['reference_n']:
            continue
        U, reference = stored[system, n], stored[system, result['reference_n']]
        delta = abs(U-reference)
        safe = abs(reference[:, 0]) > 1e-14
        result['records'].append(dict(system=system, n=n,
            scalar_relative_U_linf=float(np.max(delta[safe, 0]/abs(reference[safe, 0]))) if np.any(safe) else None,
            relative_scalar_points=int(np.sum(safe)),
            absolute_U_linf=np.max(delta, axis=0).tolist(),
            scalar_relative_U_infinity=float(np.max(delta[(coords[:, 0] == 1)&safe, 0]/abs(reference[(coords[:, 0] == 1)&safe, 0]))) if np.any((coords[:, 0] == 1)&safe) else None))
    a.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
