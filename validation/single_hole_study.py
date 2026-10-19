"""Paper single-puncture convergence and analytic near-extremal Liu Kerr.

Reproduce physical parameters and solution-reference norms of Figs.1/2.
Angular Legendre collocation replaces the paper's theta-Chebyshev grid.
Kerr counterparts are exact non-conformally-flat seeds, not BY solutions.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'python'), str(ROOT/'validation')]
from single_hole import BowenYorkSingleHole
from hispid import Backend, Hole
from physical import constraints, norms


def scaled_error(a, b):
    return float(np.max(abs(np.asarray(a)-b)/(1+abs(b))))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    b = Backend(str(a.library.resolve()))
    start = time.monotonic()
    report = dict(reference='https://arxiv.org/pdf/gr-qc/0404056',
                  scope=__doc__, library_sha256=b.library_sha256(),
                  source_sha256={str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest()
                      for f in (Path(__file__), ROOT/'python/single_hole.py', ROOT/'validation/physical.py')},
                  single_puncture=[], conformal_kerr=[], liu_extreme=[])
    def save():
        (a.output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    A, mu = np.meshgrid(np.linspace(0, 1, 6), np.cos(np.linspace(0, np.pi, 6)), indexing='ij')
    # Fixed Cartesian observers include both sides of the throat and generic
    # directions, with no singular puncture or coordinate-axis stencil.
    xyz = np.array([[.5, .1, .2], [1., .2, -.4], [3., 2., 4.], [.17, -.12, .1], [.08, .06, -.13], [8., -3., 2.]])
    direction = np.array([[.3, .4, .5], [-.7, .2, .3], [.1, -.8, .4], [.5, -.2, -.6]])
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    rh = (1+np.sqrt(1-.9999**2))/4
    xyz = np.concatenate((xyz, .9*rh*direction, 1.1*rh*direction))
    # Fixed observer, declared once; compare its exact-seed floor explicitly.
    step = np.minimum(.004, .008*np.linalg.norm(xyz, axis=1))
    report['observer'] = dict(xyz=xyz.tolist(), step=step.tolist(), order=6)
    save()
    for label, spin, momentum, reference_n in (('spin', .2, 0., 50), ('boost', 0., .2, 70)):
        ref = BowenYorkSingleHole(spin=spin, momentum=momentum, nr=reference_n, ntheta=reference_n)
        diag = ref.solve()
        assert diag['converged'], diag
        refV = ref.grid.evaluate_compactified(ref.V, A, mu)
        record = dict(case=label, spin=spin, momentum=momentum, reference_n=reference_n,
                      reference_diagnostics=diag, levels=[])
        for n in (10, 16, 24, 32, 40):
            s = BowenYorkSingleHole(spin=spin, momentum=momentum, nr=n, ntheta=n)
            diagnostics = s.solve()
            V = s.grid.evaluate_compactified(s.V, A, mu)
            physical = constraints(s.sample, xyz, step, order=6)
            row = dict(n=n, diagnostics=diagnostics,
                relative_U_linf=float(np.max(abs(1-V/refV))),
                absolute_U_linf=float(np.max(abs(V-refV))), physical=norms(physical))
            record['levels'].append(row)
            print(label, row, flush=True)
        record['reference_physical'] = norms(constraints(ref.sample, xyz, step, order=6))
        if label == 'boost':
            record['observed_order_first_four'] = float(-np.polyfit(np.log([x['n'] for x in record['levels'][:4]]),
                np.log([x['relative_U_linf'] for x in record['levels'][:4]]), 1)[0])
        # A nonzero, smooth initial guess must reach the same discrete state.
        independent = BowenYorkSingleHole(spin=spin, momentum=momentum, nr=reference_n, ntheta=reference_n)
        independent.V[:] = .003*(1-independent.grid.A[:, None])
        record['independent_start'] = dict(diagnostics=independent.solve(),
            solution_linf=float(np.max(abs(independent.V-ref.V))))
        np.savez_compressed(a.output/f'by-{label}-reference.npz', V=ref.V, A=ref.grid.A, mu=ref.grid.mu)
        report['single_puncture'].append(record)
        # An exact Kerr solution is the meaningful corresponding free metric:
        # no nonzero elliptic correction should be manufactured for this test.
        hole = Hole(1., spin=(0, 0, .2) if spin else (0, 0, 0),
                    velocity=(0, 0, momentum/np.sqrt(1+momentum*momentum)))
        sample = lambda x:b.seed(hole, x, seed_family='qi')
        seed_record = dict(case=label,
            interpretation='exact Kerr/boosted Schwarzschild seed; correction is zero',
            exact_floor=norms(constraints(sample, xyz, step, order=6)), levels=[])
        truth = {name:sample(xyz)[name].reshape(-1, 3, 3) for name in ('gamma', 'Kij')}
        for n in (10, 16, 24):
            represented = b.single_hole(chi=.2 if spin else 0., velocity=hole.velocity[2],
                                       nr=n, ntheta=n, seed_family='qi')
            fields = represented.sample(xyz)
            seed_record['levels'].append(dict(n=[n, n],
                solution_error={name:scaled_error(fields[name], truth[name]) for name in truth},
                physical=norms(constraints(represented.sample, xyz, step, order=6))))
        report['conformal_kerr'].append(seed_record)
        save()
    hole = Hole(1., spin=(0, 0, .9999))
    exact = b.seed(hole, xyz, seed_family='liu')
    exact_fields = {k:exact[k].reshape(-1, 3, 3) for k in ('gamma', 'Kij')}
    report['liu_exact_floor'] = norms(constraints(lambda x:b.seed(hole, x, seed_family='liu'), xyz, step, order=6))
    report['liu_exact'] = dict(mass=1., chi=.9999, area=8*np.pi*(1+np.sqrt(1-.9999**2)), J=.9999)
    last = None
    for n in (64, 96, 192):
        k = b.single_hole(chi=.9999, nr=n, ntheta=32, seed_family='liu')
        sampled = k.sample(xyz)
        physical = constraints(k.sample, xyz, step, order=6)
        row = dict(n=[n, 32], horizon_radius=k.horizon_radius,
                   solution_error={name:scaled_error(sampled[name], exact_fields[name]) for name in exact_fields},
                   physical=norms(physical), horizon=k.horizon())
        if last is not None:
            row['difference_from_previous'] = {name:scaled_error(sampled[name], last[name]) for name in exact_fields}
        last = sampled
        report['liu_extreme'].append(row)
        save()
        print('liu', row, flush=True)
    final = report['liu_extreme'][-1]
    report['liu_criteria'] = dict(scaled_field_linf=1e-10, mass_spin_absolute=1e-9,
        expansion_rms=1e-8, physical_rms='monotone, then within5x independently measured exact-seed FD6 floor (minimum H1e-9/M1e-11)')
    physical_passed = all(report['liu_extreme'][i+1]['physical'][name] < report['liu_extreme'][i]['physical'][name]
        for i in range(2) for name in ('H_rms', 'M_rms'))
    physical_passed = physical_passed and all(final['physical'][name] <= max(5*report['liu_exact_floor'][name], minimum)
        for name, minimum in (('H_rms', 1e-9), ('M_rms', 1e-11)))
    report['liu_passed'] = bool(max(final['solution_error'].values()) < 1e-10
        and abs(final['horizon']['mass']-1) < 1e-9 and abs(final['horizon']['chi']-.9999) < 1e-9
        and final['horizon']['expansion_rms'] < 1e-8 and physical_passed
        and all(report['liu_extreme'][i+1]['solution_error'][k] < report['liu_extreme'][i]['solution_error'][k]
                for i in range(2) for k in exact_fields))
    report['elapsed_seconds'] = time.monotonic()-start
    save()
    if not report['liu_passed']:
        raise RuntimeError('Liu analytic representation/horizon convergence gate failed')


if __name__ == '__main__':
    main()
