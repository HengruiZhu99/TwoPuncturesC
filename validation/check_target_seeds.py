"""Independent physical checks for the revised isolated seed targets.

The Cartesian finite-difference constraints and ADM integrals use only physical
fields. An additional level-set expansion check uses the physical metric's
first derivatives, with an independently specified exact horizon ellipsoid.
No elliptic binary acceptance is inferred from these zero-correction controls.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

from hispid import Backend, Hole
from configs import as_dict
from checkpoints import ROOT, library_sha
from physical import constraints, norms, charges, extrapolate
from native_loader import loaded_kokkos_images


def direction(x):
    x = np.asarray(x, dtype=float)
    return x / np.linalg.norm(x)


def write_progress(path,bound,metadata,output):
    """Atomic update of a reserved new attempt; changed inputs leave old evidence."""
    if any(hashlib.sha256(Path(file).read_bytes()).hexdigest()!=sha for file,sha in bound.items()):
        raise ValueError('bound exact-seed control implementation or native image changed')
    payload=dict(output,bound_artifacts_sha256=bound,**metadata)
    temporary=path.with_name(path.name+'.tmp')
    with temporary.open('x') as stream:stream.write(json.dumps(payload,indent=2)+'\n')
    os.replace(temporary,path)


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


def run(backend, radii, targets=None, quadrature_levels=None, align_polar_axis=False,progress=None,seed_execution='reference',step_factors=(2,1,.5),all_angular_changes=False):
    axes = np.array([[.73, .31, .61], [-.41, .82, .39], [.22, -.51, .83],
                     [-.69, -.44, .57], [.39, .73, -.56], [.81, -.38, -.45]])
    axes /= np.linalg.norm(axes, axis=1)[:, None]
    spin_axis, boost_axis = direction([.2, -.3, .4]), direction([-.7, .2, .3])
    cases = targets or [('spin95', .95, 0), ('boost885', 0, .885),
                       ('spin95_boost885_generic', .95, .885)]
    quadrature_levels=quadrature_levels or ((16,32),(24,48),(32,64))
    output = dict(library_sha256=library_sha(backend), cases=[], passed=False,completed=False,
                  seed_execution=seed_execution,exact_horizon_geometry='host',
                  criteria=dict(physical_constraint_rms=1e-7,
                                ADM_absolute_error=1e-5,
                                ADM_quadrature_change=1e-5,
                                ADM_radial_fit_change=1e-5,
                                exact_horizon_expansion_max=1e-10))
    def save():
        if library_sha(backend)!=output['library_sha256']:
            raise ValueError('exact-seed native image changed')
        if progress:progress(output)
    save()
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
        frame=None
        if align_polar_axis:
            polar=direction(v if speed else S)
            transverse=np.eye(3)[np.argmin(abs(polar))]
            transverse=direction(transverse-(transverse@polar)*polar)
            frame=np.column_stack([polar,transverse,np.cross(polar,transverse)])
        rh = .5 * np.sqrt(1 - chi**2)
        ray_radius = rh / np.sqrt(1 + G**2 * (axes @ v)**2)
        near = np.concatenate([factor * ray_radius[:, None] * axes
                               for factor in (1.5, 3)] + [.6 * axes, axes])
        xyz = np.r_[near, 4 * axes, 8 * axes, 16 * axes]
        step = np.minimum(.002, .001 * np.linalg.norm(xyz, axis=1))
        sample = lambda x: backend.seed(hole, x, config.conformal_choice,execution=seed_execution)
        sequence = []
        expected = np.r_[G, G * v, G * S - G**2 / (G + 1) * (v @ S) * v]
        quadratures = []
        record = dict(case=label, config=as_dict(config), seed_rest_chi=chi, lab_speed=speed,
                      input_lorentz_factor=float(G),polar_frame=frame.tolist() if frame is not None else None,
                      xyz=xyz.tolist(),verifier_steps=step.tolist(),constraint_sequence=sequence,
                      radii=list(radii),charge_quadratures=quadratures,expected_charges=expected.tolist(),
                      passed=False,completed=False)
        output['cases'].append(record);save()
        for factor in step_factors:
            residual = constraints(sample, xyz, factor * step)
            sequence.append(dict(step_factor=factor,
                                 near=norms(residual, np.arange(len(xyz)) < len(near)),
                                 bulk=norms(residual, np.arange(len(xyz)) >= len(near))))
            save()
        with backend.create_sampler(config) as solution:
            expansion = horizon_expansion(solution, hole, axes)
        record['exact_horizon_expansion']=expansion.tolist();save()
        for nt, np_ in quadrature_levels:
            q=[];quadrature=dict(ntheta=nt,nphi=np_,charges=[],radial_fits=[],completed=False)
            quadratures.append(quadrature);save()
            for r in radii:
                value=charges(sample,r,ntheta=nt,nphi=np_,polar_frame=frame)
                q.append(value);quadrature['charges'].append(np.asarray(value).tolist());save()
            fits = [extrapolate(radii[i:i+4], q[i:i+4]) for i in range(len(radii)-3)]
            quadrature.update(radial_fits=np.array(fits).tolist(),completed=True);save()
        final = np.array(quadratures[-1]['radial_fits'][-1])
        error = np.max(abs(final - expected))
        quad_change = np.max(abs(final - quadratures[-2]['radial_fits'][-1]))
        adjacent_changes=[float(np.max(abs(np.asarray(b['radial_fits'][-1])-a['radial_fits'][-1])))
                          for a,b in zip(quadratures,quadratures[1:])]
        radial_change = np.max(abs(final - quadratures[-1]['radial_fits'][-2]))
        constraint_pass = all(sequence[-1][region][component] < 1e-7
                              for region in ('near', 'bulk') for component in ('H_rms', 'M_rms'))
        passed = (constraint_pass and error < 1e-5 and quad_change < 1e-5
                  and radial_change < 1e-5 and np.max(abs(expansion)) < 1e-10)
        if all_angular_changes:passed &= max(adjacent_changes)<1e-5
        record.update(charge_error_max=float(error), quadrature_change_max=float(quad_change),
                      adjacent_quadrature_change_max=adjacent_changes,all_angular_changes_required=all_angular_changes,
                      radial_fit_change_max=float(radial_change), passed=bool(passed),
                      seconds=time.monotonic() - start,completed=True)
        save()
        print(label, sequence[-1], 'charge error', error, 'Theta', expansion, flush=True)
    output['passed'] = all(case['passed'] for case in output['cases'])
    output['completed']=True;save()
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--library', required=True)
    parser.add_argument('--output', default='validation/target_seed_controls.json')
    parser.add_argument('--radii', default='40,80,160,320,640,1280,2560,5120,10240')
    parser.add_argument('--extreme',action='store_true',help='fresh separate chi=.99/Gamma=10 controls, with refined beam-aligned charge quadrature')
    parser.add_argument('--target-case',choices=('spin99','gamma10'),help='one separately retained extreme seed refinement; requires --extreme')
    parser.add_argument('--quadratures',help='explicit ntheta:nphi levels; defaults remain unchanged')
    parser.add_argument('--step-factors',default='2,1,.5',help='strictly decreasing positive Cartesian verifier step factors')
    parser.add_argument('--seed-execution',choices=('reference','kokkos'),default='reference',help='evaluator used by independent Cartesian FD and charge quadrature; exact-horizon gradient check remains host')
    parser.add_argument('--threads',type=int,default=1)
    args = parser.parse_args()
    if args.target_case and not args.extreme:raise ValueError('target-specific refinement requires --extreme')
    step_factors=list(map(float,args.step_factors.split(',')))
    if (len(step_factors)<3 or not np.isfinite(step_factors).all() or min(step_factors)<=0
        or not np.all(np.diff(step_factors)<0)):
        raise ValueError('at least three finite positive decreasing verifier steps required')
    quadratures=[tuple(map(int,x.split(':'))) for x in args.quadratures.split(',')] if args.quadratures else None
    if quadratures is not None and (len(quadratures)<2 or any(len(x)!=2 or min(x)<8 for x in quadratures)
        or len(set(quadratures))!=len(quadratures)
        or any(b[0]<a[0] or b[1]<a[1] for a,b in zip(quadratures,quadratures[1:]))):
        raise ValueError('distinct nondecreasing angular quadratures of at least8 points required')
    if not 1<=args.threads<=16:raise ValueError('between1 and16 host threads required')
    radii = list(map(float, args.radii.split(',')))
    if len(radii)<5 or min(radii)<=0 or not np.isfinite(radii).all() or not np.all(np.diff(radii)>0):
        raise ValueError('at least five finite positive increasing charge radii are required')
    path=(ROOT/args.output).resolve();path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as stream:
        stream.write(json.dumps(dict(passed=False,completed=False,stage='loading_bound_library'))+'\n')
    latest={};bound=None;metadata=None
    def progress(output):
        latest['output']=output
        write_progress(path,bound,metadata,output)
    try:
        backend=Backend(args.library)
        from execution import select,name,concurrency,device_description
        select(backend.lib,args.seed_execution,args.threads)
        device=device_description(backend.lib) if args.seed_execution=='kokkos' else None
        if args.seed_execution=='kokkos' and name(backend.lib)=='Cuda' and (not device or device['visible_count']!=1 or device['visible_ordinal']!=0):
            raise ValueError('execution seed controls require one visible CUDA device')
        paths=[Path(__file__).resolve()]+[Path(__import__(name).__file__).resolve() for name in
            ('hispid','configs','checkpoints','physical','native_loader','execution')]
        bound={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
        bound[str(backend.path)]=library_sha(backend);bound.update(backend.dependency_images)
        runtime_images=loaded_kokkos_images();bound.update(runtime_images)
        metadata=dict(library_dependency_images=backend.dependency_images,
            control_scope='requested_target' if args.target_case else 'full',target_case=args.target_case,
            verifier_step_factors=step_factors,all_angular_changes_required=bool(args.target_case),
            unknown_parameterization_id=backend.parameterization(),collocation_maps=backend.parameterization_maps(),
            seed_execution=args.seed_execution,compiled_execution=name(backend.lib) if args.seed_execution=='kokkos' else 'reference',
            execution_concurrency=concurrency(backend.lib) if args.seed_execution=='kokkos' else 1,device=device,
            seed_scalar_digits=53 if args.seed_execution=='kokkos' else None,exact_horizon_geometry='host',
            runtime_images=runtime_images,native_images={str(backend.path):library_sha(backend),**backend.dependency_images,**runtime_images})
        targets=[('spin99',.99,0),('gamma10',0,float(np.sqrt(.99)))] if args.extreme else None
        if args.target_case:targets=[case for case in targets if case[0]==args.target_case]
        result = run(backend, radii,targets=targets,
            quadrature_levels=quadratures or (((64,64),(128,128),(192,192)) if args.extreme else None),
            align_polar_axis=args.extreme,progress=progress,seed_execution=args.seed_execution,step_factors=step_factors,
            all_angular_changes=bool(args.target_case))
    except BaseException as error:
        output=latest.get('output',dict(stage='setup_failed',cases=[]))
        output.update(passed=False,completed=False,failure=dict(type=type(error).__name__,message=str(error)))
        try:write_progress(path,bound or {},metadata or {},output)
        except BaseException as write_error:
            print('Failure record unavailable; prior snapshot preserved: '+str(write_error),file=sys.stderr)
        raise
    raise SystemExit(0 if result['passed'] else 1)
