"""Distinct asymptotic seed gate: independent raw-tensor ADM surface integrals.

One host seed path; no elliptic solve or backend timing matrix. Charge observer
uses Cartesian finite differences, not the seed's analytic metric derivatives.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'python'))
from hispid import Backend, Hole
from physical import charges, extrapolate


def run(library, output, only=None, radius_scale=1., reuse=None):
    backend = Backend(str(library.resolve()))
    if not np.isfinite(radius_scale) or radius_scale <= 0: raise ValueError('positive finite radius scale required')
    radii = [radius_scale*r for r in (256., 512., 1024., 2048.)]
    result = dict(kind='trumpet_independent_ADM', library_sha256=backend.loaded_sha256,
                  radii=radii, relative_tolerance=1e-4, cases=[], completed=False,
                  binary_acceptance=False)
    root = Path(__file__).resolve().parents[1]
    result['source_sha256'] = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [Path(__file__).resolve(), root/'validation/physical.py',
                  root/'python/hispid.py', root/'src/HiSpID_trumpet.hpp']}
    cached = json.loads(reuse.read_text()) if reuse else None
    if cached:
        assert cached["library_sha256"] == backend.loaded_sha256
        assert cached["source_sha256"]["validation/physical.py"] == result["source_sha256"]["validation/physical.py"]
        result["reused_result"] = dict(path=str(reuse.resolve()), sha256=hashlib.sha256(reuse.read_bytes()).hexdigest())
    def save():
        output.write_text(json.dumps(result, indent=2)+'\n')
    for name, spin, velocity, levels in [
        ('spin99', (0, 0, .99), (0, 0, 0), [(16, 32), (24, 48)]),
        ('generic', (.2, -.3, .7), (.25, .1, -.2), [(16, 32), (24, 48)]),
        ('gamma10', (0, 0, 0), (np.sqrt(.99), 0, 0), [(96, 8), (128, 8)])]:
        if only and name != only: continue
        start = time.monotonic()
        hole = Hole(1, spin=spin, velocity=velocity)
        v, S = np.array(velocity), np.array(spin)
        G = 1/np.sqrt(1-v@v)
        expected = np.r_[G, G*v, G*S-G**2/(G+1)*(v@S)*v]
        axis = v if np.linalg.norm(v) else S
        axis = axis/np.linalg.norm(axis)
        transverse = np.eye(3)[np.argmin(abs(axis))]
        transverse -= (transverse@axis)*axis
        transverse /= np.linalg.norm(transverse)
        frame = np.column_stack([axis, transverse, np.cross(axis, transverse)])
        row = dict(case=name, spin=spin, velocity=list(velocity), expected=expected.tolist(),
                   quadratures=[], completed=False)
        result['cases'].append(row); save()
        sample = lambda x: backend.seed(hole, x, seed_family='trumpet_r0_m')
        for nt, nf in levels:
            q = dict(ntheta=nt, nphi=nf, values=[])
            row['quadratures'].append(q)
            for r in radii:
                prior = None
                if cached and r in cached['radii']:
                    for cr in cached['cases']:
                        if cr['case'] == name and cr['spin'] == list(spin) and cr['velocity'] == list(velocity):
                            for cq in cr['quadratures']:
                                if cq['ntheta'] == nt and cq['nphi'] == nf:
                                    prior = cq['values'][cached['radii'].index(r)]
                value = prior if prior is not None else charges(sample, r, ntheta=nt, nphi=nf, polar_frame=frame).tolist()
                if not np.all(np.isfinite(value)): raise ValueError('nonfinite ADM integral')
                q['values'].append(value); save()
            q['fits'] = [extrapolate(radii[i:i+3], q['values'][i:i+3]).tolist()
                         for i in (0, 1)]
            save()
        final = np.array(row['quadratures'][-1]['fits'][-1])
        # One mass unit is the absolute scale for a vanishing charge component.
        scale = np.maximum(1., abs(expected))
        row['scaled_error'] = float(np.max(abs(final-expected)/scale))
        row['angular_change'] = float(np.max(abs(final-row['quadratures'][0]['fits'][-1])/scale))
        row['radial_change'] = float(np.max(abs(final-row['quadratures'][-1]['fits'][0])/scale))
        row['passed'] = max(row[k] for k in ('scaled_error','angular_change','radial_change')) < 1e-4
        row['seconds'] = time.monotonic()-start
        row['completed'] = True; save()
        print(name, row['scaled_error'], row['angular_change'], row['radial_change'], flush=True)
    result['completed'] = True
    result['passed'] = all(row['passed'] for row in result['cases'])
    backend.library_sha256(); save()
    return result

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--library', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--case', choices=['spin99','generic','gamma10'])
    ap.add_argument('--radius-scale', type=float, default=1.)
    ap.add_argument('--reuse', type=Path)
    args = ap.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    run(args.library, args.output, args.case, args.radius_scale, args.reuse)
