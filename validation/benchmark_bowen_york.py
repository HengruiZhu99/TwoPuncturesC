"""Sequential cold-start BY/HiSpID costs at matched counts and bare inputs.

This does not equate seed geometries, horizon properties or residual norms.
"""
import argparse
import ctypes as C
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
from configs import as_dict, moderate, target_binary
from hispid import Backend

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_image(lib, symbol, expected):
    if sys.platform != 'darwin' and not sys.platform.startswith('linux'):
        raise RuntimeError('loaded image verification requires dladdr')
    class DlInfo(C.Structure):
        _fields_ = [('filename', C.c_char_p), ('base', C.c_void_p),
                    ('symbol', C.c_char_p), ('symbol_address', C.c_void_p)]
    dladdr = C.CDLL(None).dladdr
    dladdr.argtypes = [C.c_void_p, C.POINTER(DlInfo)]
    dladdr.restype = C.c_int
    info = DlInfo()
    address = C.cast(getattr(lib, symbol), C.c_void_p).value
    if not dladdr(address, C.byref(info)) or not info.filename:
        raise RuntimeError('cannot identify loaded benchmark library')
    if Path(info.filename.decode()).resolve(strict=True) != expected:
        raise RuntimeError('benchmark native loader reused a different image')


def configuration(backend, case, grid, tolerance):
    if case == 'moderate':
        config = moderate(backend, grid[0], grid[2]); config.far_radius = 0
    else:
        config = target_binary(backend, grid[0], grid[2],
                               speed=.885 if case == 'boost885' else 0,
                               spin=.95 if case == 'spin95' else 0)
    config.n[:] = grid
    config.tolerance = tolerance
    config.max_newton = 24
    config.max_krylov = 2000
    config.krylov_restart = 64
    config.memory_limit_mib = 8192
    return config


def worker(args):
    # Only HiSpID workers load the HiSpID library. BY consumes a saved input.
    if args.mode == 'hispid':
        backend = Backend(str(Path(args.hispid_library).resolve(strict=True)))
        config = configuration(backend, args.case, args.grid, args.tolerance)
        start = time.monotonic()
        with backend.create(config) as data:
            created = time.monotonic()
            diagnostics = data.solve()
            finished = time.monotonic()
            data.sample([[0., 2., 1.]])
            ready = time.monotonic()
        result = dict(config=as_dict(config), diagnostics=diagnostics,
                      creation_seconds=created-start,
                      setup_and_solve_seconds=finished-start,
                      solve_seconds=diagnostics['seconds'],
                      first_sample_seconds=ready-finished,
                      ready_to_sample_seconds=ready-start,
                      converged=diagnostics['status'] == 0,
                      library_sha256=backend.library_sha256(),
                      loaded_image_verified=True)
    else:
        config = json.loads(Path(args.input).read_text())
        path = Path(args.by_library).resolve(strict=True)
        expected_sha = digest(path)
        lib = C.CDLL(str(path))
        verify_image(lib, 'TwoPunctures_make_initial_data', path)
        verify_image(lib, 'benchmark_newton_seconds', path)
        signatures = {
            'TwoPunctures_params_set_default': (None, []),
            'TwoPunctures_params_set_Real': (None, [C.c_char_p, C.c_double]),
            'TwoPunctures_params_set_Int': (None, [C.c_char_p, C.c_int]),
            'TwoPunctures_make_initial_data': (C.c_void_p, []),
            'TwoPunctures_diagnostics': (C.c_int, [C.c_void_p,
                C.POINTER(C.c_double), C.POINTER(C.c_double), C.POINTER(C.c_double)]),
            'TwoPunctures_finalise': (None, [C.c_void_p]),
            'TwoPunctures_sample_points': (C.c_int, [C.c_void_p, C.c_int,
                C.POINTER(C.c_double), C.POINTER(C.c_double), C.POINTER(C.c_double),
                C.POINTER(C.c_double), C.POINTER(C.c_double)]),
            'benchmark_newton_seconds': (C.c_double, []),
        }
        for name, (ret, params) in signatures.items():
            getattr(lib, name).restype = ret
            getattr(lib, name).argtypes = params
        lib.TwoPunctures_params_set_default()
        reals = {'par_b': config['hole'][0]['center'][0],
                 'par_m_plus': config['hole'][0]['mass'],
                 'par_m_minus': config['hole'][1]['mass'],
                 'Newton_tol': config['tolerance']}
        momenta, spins = [], []
        for side, hole in zip(('plus', 'minus'), config['hole']):
            velocity = np.array(hole['velocity'])
            lorentz = 1/np.sqrt(1-velocity@velocity)
            momentum = hole['mass'] * lorentz * velocity
            rest_spin = np.array(hole['spin'])
            spin = lorentz*rest_spin - lorentz**2/(1+lorentz)*velocity*(velocity@rest_spin)
            momenta.append(momentum.tolist())
            spins.append(spin.tolist())
            for axis in range(3):
                reals[f'par_P_{side}{axis+1}'] = float(momentum[axis])
                reals[f'par_S_{side}{axis+1}'] = float(spin[axis])
        integers = dict(zip(('npoints_A', 'npoints_B', 'npoints_phi'), config['n']))
        integers.update(Newton_maxit=config['max_newton'], give_bare_mass=1,
                        use_external_initial_guess=0, solve_momentum_constraint=0,
                        grid_setup_method=1,
                        verbose=int(args.by_verbose))
        for name, value in reals.items():
            lib.TwoPunctures_params_set_Real(name.encode(), value)
        for name, value in integers.items():
            lib.TwoPunctures_params_set_Int(name.encode(), value)
        start = time.monotonic()
        data = lib.TwoPunctures_make_initial_data()
        finished = time.monotonic()
        if not data:
            raise RuntimeError('BY allocation/solve failed')
        xyz = (C.c_double*3)(0., 2., 1.)
        lapse, psi = (C.c_double*1)(), (C.c_double*1)()
        gamma, curvature = (C.c_double*6)(), (C.c_double*6)()
        sample_status = lib.TwoPunctures_sample_points(data, 1, xyz, lapse, psi, gamma, curvature)
        ready = time.monotonic()
        if sample_status:
            raise RuntimeError('BY first sampling failed')
        residual, energy = C.c_double(), C.c_double()
        masses = (C.c_double*2)()
        status = lib.TwoPunctures_diagnostics(data, C.byref(residual), C.byref(energy), masses)
        lib.TwoPunctures_finalise(data)
        if digest(path) != expected_sha:
            raise RuntimeError('BY library changed on disk during benchmark')
        result = dict(config=config, by_real_parameters=reals,
                      by_integer_parameters=integers, by_momenta=momenta,
                      by_lab_intrinsic_spins=spins,
                      solve_seconds=lib.benchmark_newton_seconds(),
                      setup_and_solve_seconds=finished-start,
                      first_sample_seconds=ready-finished,
                      ready_to_sample_seconds=ready-start,
                      diagnostics=dict(status=status, residual_linf=residual.value,
                                       adm_energy=energy.value, puncture_masses=list(masses)),
                      converged=bool(status == 0 and np.isfinite(residual.value)
                                     and residual.value <= config['tolerance']),
                      library_sha256=expected_sha, loaded_image_verified=True,
                      sample_method='spectral',
                      verbose_inside_timer=args.by_verbose)
    result.update(case=args.case, cpu_threads=1, initial_guess='zero',
                  max_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                  * (1 if sys.platform == 'darwin' else 1024))
    if args.mode == 'hispid':
        result['sample_method'] = 'spectral'
    Path(args.worker_output).write_text(json.dumps(result, indent=2)+'\n')
    print('BENCHMARK', json.dumps(result), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hispid-library', required=True)
    parser.add_argument('--by-library', required=True)
    parser.add_argument('--case', choices=('spin95', 'boost885', 'moderate'), default='spin95')
    parser.add_argument('--grid', type=lambda x: tuple(map(int, x.split(':'))), default=(80,160,16))
    parser.add_argument('--tolerance', type=float, default=1e-10)
    parser.add_argument('--output', required=True)
    parser.add_argument('--timeout', type=float, default=1200)
    parser.add_argument('--mode', choices=('by', 'hispid'))
    parser.add_argument('--input')
    parser.add_argument('--worker-output')
    parser.add_argument('--by-verbose', action='store_true')
    args = parser.parse_args()
    if len(args.grid) != 3 or min(args.grid) < 4 or args.grid[2] % 2:
        parser.error('grid needs three sizes >=4 and an even Fourier size')
    if not np.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error('positive finite tolerance required')
    if args.mode:
        worker(args); return
    output = Path(args.output)
    raw = ROOT/'validation/raw'/output.stem
    if output.exists():
        raise FileExistsError(output)
    raw.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    records = {}
    for mode in ('hispid', 'by'):
        path = raw/(mode+'.json')
        cmd = [sys.executable, str(Path(__file__).resolve()),
               '--hispid-library', args.hispid_library, '--by-library', args.by_library,
               '--case', args.case, '--grid', ':'.join(map(str, args.grid)),
               '--tolerance', str(args.tolerance), '--output', str(output),
               '--mode', mode, '--worker-output', str(path)]
        if args.by_verbose:
            cmd += ['--by-verbose']
        if mode == 'by':
            input_path = raw/'input.json'
            input_path.write_text(json.dumps(records['hispid']['config'], indent=2)+'\n')
            cmd += ['--input', str(input_path)]
        print(f'Starting {mode} {args.case} {args.grid}', flush=True)
        with (raw/(mode+'.log')).open('w') as log:
            subprocess.run(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT,
                           timeout=args.timeout, check=True)
        records[mode] = json.loads(path.read_text())
        print(mode, records[mode]['solve_seconds'], records[mode]['converged'], flush=True)
    converged = all(r['converged'] for r in records.values())
    result = dict(case=args.case, grid=args.grid, records=records,
                  both_native_stopping_criteria_met=converged,
                  cold_start=True, cpu_threads=1, isolated_sequential_processes=True,
                  hispid_over_by_solve_ratio=records['hispid']['solve_seconds']/records['by']['solve_seconds'] if converged else None,
                  hispid_over_by_setup_and_solve_ratio=records['hispid']['setup_and_solve_seconds']/records['by']['setup_and_solve_seconds'] if converged else None,
                  hispid_over_by_ready_to_sample_ratio=records['hispid']['ready_to_sample_seconds']/records['by']['ready_to_sample_seconds'] if converged else None,
                  note='Matched collocation counts, bare masses and centers; BY P=m*Gamma*v and S=Gamma*Srest-Gamma^2/(1+Gamma)*v*(v dot Srest) match the isolated seed lab-frame intrinsic charges. Native stopping tolerances have the same numeric value on each native weighted equation residual. BY solves one Hamiltonian equation with analytic momentum; HiSpID solves four coupled curved equations. Maps/seed geometries, physical accuracy, final horizon masses/spins and velocities are not equated. One paired run is not a scaling curve or binary validation.',
                  validation_acceptance_inherited=False)
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
