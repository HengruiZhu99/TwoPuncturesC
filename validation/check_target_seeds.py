"""Independent physical checks for the revised isolated seed targets.

The Cartesian finite-difference constraints and ADM integrals use only physical
fields. An additional level-set expansion check uses the physical metric's
first derivatives, with an independently specified exact horizon ellipsoid.
No elliptic binary acceptance is inferred from these zero-correction controls.
"""
import argparse
import json
import time

import numpy as np

from hispid import Backend, Hole
from configs import as_dict
from checkpoints import ROOT, library_sha
from physical import constraints, norms, charges, extrapolate


def direction(x):
    x = np.asarray(x, dtype=float)
    return x / np.linalg.norm(x)


def horizon_expansion(solution, hole, directions):
    """Theta of F=X^T Q X-r_h^2, including its covariant Hessian."""
    v = np.array(hole.velocity)
    chi = np.linalg.norm(hole.spin) / hole.mass**2
    G = 1 / np.sqrt(1 - v @ v)
    Q = np.eye(3) + G**2 * np.outer(v, v)
    rh = .5 * hole.mass * np.sqrt(1 - chi**2)
    radius = rh / np.sqrt(np.einsum('ni,ij,nj->n', directions, Q, directions))
    X = radius[:, None] * directions
    data = solution.sample_with_derivatives(X + np.array(hole.center))
    metric = data['gamma'].reshape(-1, 3, 3)
    K, dg = data['Kij'].reshape(-1, 3, 3), data['dgamma']
    inverse = np.linalg.inv(metric)
    gradient = 2 * X @ Q
    length = np.sqrt(np.einsum('ni,nij,nj->n', gradient, inverse, gradient))
    normal = np.einsum('nij,nj->ni', inverse, gradient) / length[:, None]
    connection = np.zeros((len(X), 3, 3, 3))
    for k in range(3):
        for i in range(3):
            for j in range(3):
                connection[:, k, i, j] = .5 * np.einsum(
                    'nl,nl->n', inverse[:, k],
                    dg[:, i, j] + dg[:, j, i] - dg[:, :, i, j])
    hessian = 2 * Q - np.einsum('nkij,nk->nij', connection, gradient)
    projector = inverse - np.einsum('ni,nj->nij', normal, normal)
    return (np.einsum('nij,nij->n', projector, hessian) / length
            + np.einsum('ni,nij,nj->n', normal, K, normal)
            - np.einsum('nij,nij->n', inverse, K))


def run(backend, radii):
    axes = np.array([[.73, .31, .61], [-.41, .82, .39], [.22, -.51, .83],
                     [-.69, -.44, .57], [.39, .73, -.56], [.81, -.38, -.45]])
    axes /= np.linalg.norm(axes, axis=1)[:, None]
    spin_axis, boost_axis = direction([.2, -.3, .4]), direction([-.7, .2, .3])
    cases = [('spin95', .95, 0), ('boost885', 0, .885),
             ('spin95_boost885_generic', .95, .885)]
    output = dict(library_sha256=library_sha(backend), cases=[], passed=False,
                  criteria=dict(physical_constraint_rms=1e-7,
                                ADM_absolute_error=1e-5,
                                ADM_quadrature_change=1e-5,
                                ADM_radial_fit_change=1e-5,
                                exact_horizon_expansion_max=1e-10))
    for label, chi, speed in cases:
        start = time.monotonic()
        hole = Hole(1, spin=chi * spin_axis, velocity=speed * boost_axis)
        config = backend.config()
        config.n[:] = [6, 6, 4]
        config.hole[0] = hole
        config.hole[1] = Hole(0, center=(-6, 0, 0))
        config.omega[:] = [0, 0]
        config.far_radius = 0
        config.inner_min[:] = [0, 0]
        config.inner_max[:] = [0, 0]
        config.inner_flatten = 0
        config.conformal_choice = 0
        v, S = np.array(hole.velocity), np.array(hole.spin)
        G = 1 / np.sqrt(1 - v @ v)
        rh = .5 * np.sqrt(1 - chi**2)
        ray_radius = rh / np.sqrt(1 + G**2 * (axes @ v)**2)
        near = np.concatenate([factor * ray_radius[:, None] * axes
                               for factor in (1.5, 3)] + [.6 * axes, axes])
        xyz = np.r_[near, 4 * axes, 8 * axes, 16 * axes]
        step = np.minimum(.002, .001 * np.linalg.norm(xyz, axis=1))
        sample = lambda x: backend.seed(hole, x, config.conformal_choice)
        sequence = []
        for factor in (2, 1, .5):
            residual = constraints(sample, xyz, factor * step)
            sequence.append(dict(step_factor=factor,
                                 near=norms(residual, np.arange(len(xyz)) < len(near)),
                                 bulk=norms(residual, np.arange(len(xyz)) >= len(near))))
        with backend.create_sampler(config) as solution:
            expansion = horizon_expansion(solution, hole, axes)
        expected = np.r_[G, G * v, G * S - G**2 / (G + 1) * (v @ S) * v]
        quadratures = []
        for nt, np_ in ((16, 32), (24, 48), (32, 64)):
            q = [charges(sample, r, ntheta=nt, nphi=np_) for r in radii]
            fits = [extrapolate(radii[i:i+4], q[i:i+4]) for i in range(len(radii)-3)]
            quadratures.append(dict(ntheta=nt, nphi=np_, charges=np.array(q).tolist(),
                                    radial_fits=np.array(fits).tolist()))
        final = np.array(quadratures[-1]['radial_fits'][-1])
        error = np.max(abs(final - expected))
        quad_change = np.max(abs(final - quadratures[-2]['radial_fits'][-1]))
        radial_change = np.max(abs(final - quadratures[-1]['radial_fits'][-2]))
        constraint_pass = all(sequence[-1][region][component] < 1e-7
                              for region in ('near', 'bulk') for component in ('H_rms', 'M_rms'))
        passed = (constraint_pass and error < 1e-5 and quad_change < 1e-5
                  and radial_change < 1e-5 and np.max(abs(expansion)) < 1e-10)
        record = dict(case=label, config=as_dict(config), seed_rest_chi=chi, lab_speed=speed,
                      xyz=xyz.tolist(), verifier_steps=step.tolist(), constraint_sequence=sequence,
                      exact_horizon_expansion=expansion.tolist(), radii=radii,
                      charge_quadratures=quadratures, expected_charges=expected.tolist(),
                      charge_error_max=float(error), quadrature_change_max=float(quad_change),
                      radial_fit_change_max=float(radial_change), passed=bool(passed),
                      seconds=time.monotonic() - start)
        output['cases'].append(record)
        print(label, sequence[-1], 'charge error', error, 'Theta', expansion, flush=True)
    output['passed'] = all(case['passed'] for case in output['cases'])
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--library', required=True)
    parser.add_argument('--output', default='validation/target_seed_controls.json')
    parser.add_argument('--radii', default='40,80,160,320,640,1280,2560,5120,10240')
    args = parser.parse_args()
    radii = list(map(float, args.radii.split(',')))
    if len(radii)<5 or min(radii)<=0 or not np.isfinite(radii).all() or not np.all(np.diff(radii)>0):
        raise ValueError('at least five finite positive increasing charge radii are required')
    result = run(Backend(args.library), radii)
    (ROOT / args.output).write_text(json.dumps(result, indent=2) + '\n')
    raise SystemExit(0 if result['passed'] else 1)
