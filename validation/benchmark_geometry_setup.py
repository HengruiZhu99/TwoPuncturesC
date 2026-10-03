"""Supplemental, paired HiSpID setup costs in an unchanged native image.

This measures construction and first sampling, without a nonlinear solve.
Run sequentially after the frozen solve campaign, in a one-GPU Slurm step.
The standalone worker also permits small local controls; those are not the
declared Perlmutter performance matrix or physical-binary qualification.
"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / name) for name in ('python', 'examples', 'validation')]
from benchmark_bowen_york import SAMPLE_POINTS, configuration, digest
from benchmark_kokkos import GPUProcessMemory, hardware, hardware_class, manifest_images, verify_artifacts
from configs import as_dict
from execution import concurrency, device_description, name, reset_statistics, select, statistics
from hispid import Backend, Config, Hole, Point
from native_loader import loaded_kokkos_images

GRIDS = [[40, 80, 16], [80, 160, 16], [128, 256, 28]]
VARIANTS = {'serial', 'openmp1', 'openmp2', 'openmp4', 'openmp8', 'openmp16', 'cuda'}
SETUP_TESTS = {'hispid_kokkos_setup', 'hispid_execution_seed_export'}
BOUND = 1e-10


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def frozen_arrays(grid):
    """Small deterministic modal unknowns and a distinct JVP direction."""
    na, nb, nphi = grid
    a = np.cos(np.pi * (np.arange(na) + .5) / na)[None, None, :]
    b = np.cos(np.pi * (np.arange(nb) + .5) / nb)[None, :, None]
    phi = (2 * np.pi * np.arange(nphi) / nphi)[:, None, None]
    envelope = (1 - a*a)**2 * (1 - b*b)**2
    fields = [envelope * (1 + .2*a + .1*b) * np.cos((q % 3 + 1)*phi)
              for q in range(4)]
    direction = np.stack(fields, axis=-1).reshape(-1) * 1e-3
    unknowns = np.stack([np.broadcast_to(envelope * (1 + .1*a*b) *
                         (1 + .1*np.sin((q % 3 + 1)*phi)), (nphi, nb, na))
                         for q in range(4)], axis=-1).reshape(-1) * 1e-5
    return dict(direction=direction, unknowns=unknowns, points=SAMPLE_POINTS.copy())


def decode_config(values):
    config = Config()
    if set(values) != {name for name, _ in Config._fields_}:
        raise ValueError('frozen configuration has unexpected fields')
    for key, value in values.items():
        if key == 'hole':
            for i, hole in enumerate(value):
                config.hole[i] = Hole(**hole)
        elif isinstance(value, list):
            getattr(config, key)[:] = value
        else:
            setattr(config, key, value)
    if as_dict(config) != values:
        raise ValueError('frozen configuration did not decode exactly')
    return config


def rss_bytes():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(rss if sys.platform == 'darwin' else rss * 1024)


def worker(args):
    # Input decoding and image loading precede runtime initialization/timers.
    config_bytes, array_bytes = Path(args.input).read_bytes(), Path(args.input_arrays).read_bytes()
    config_hash, arrays_hash = hashlib.sha256(config_bytes).hexdigest(), hashlib.sha256(array_bytes).hexdigest()
    values = json.loads(config_bytes)
    with np.load(io.BytesIO(array_bytes), allow_pickle=False) as saved:
        inputs = {key: saved[key].copy() for key in ('direction', 'unknowns', 'points')}
    del array_bytes, config_bytes
    config = decode_config(values)
    size = 4 * int(np.prod(config.n))
    if any(inputs[key].shape != (size,) for key in ('direction', 'unknowns')):
        raise ValueError('frozen arrays have incorrect unknown counts')
    if inputs['points'].shape != SAMPLE_POINTS.shape or not all(np.all(np.isfinite(v)) for v in inputs.values()):
        raise ValueError('nonfinite or incorrectly shaped frozen inputs')
    backend = Backend(str(Path(args.library).resolve(strict=True)))
    start = time.monotonic()
    select(backend.lib, 'kokkos', args.threads)
    initialization_seconds = time.monotonic() - start
    initialization_statistics = statistics(backend.lib)
    initialization_rss = rss_bytes()
    reset_statistics(backend.lib)
    actual_space = name(backend.lib)
    actual_threads = concurrency(backend.lib)
    device = device_description(backend.lib)
    start = time.monotonic()
    with backend.create(config, execution='kokkos', geometry=args.geometry) as data:
        created = time.monotonic()
        constructor_statistics = statistics(backend.lib)
        constructor_setup = data.setup_statistics()
        constructor_rss = rss_bytes()
        sampled_at = time.monotonic()
        first_sample = data.sample(inputs['points'])
        sampled = time.monotonic()
        # Capture the timed phase before allocating verification/state arrays.
        measured_rss = rss_bytes()
        setup = data.setup_statistics()
        measured_statistics = statistics(backend.lib)
        snapshot = dict(initialization_seconds=initialization_seconds,
                        creation_seconds=created-start,
                        first_sample_seconds=sampled-sampled_at,
                        ready_to_sample_seconds=(created-start) + (sampled-sampled_at),
                        wall_creation_through_sample_seconds=sampled-start,
                        instrumentation_gap_seconds=sampled_at-created,
                        initialization_rss_bytes=initialization_rss,
                        constructor_rss_bytes=constructor_rss,
                        max_rss_bytes=measured_rss,
                        initialization_statistics=initialization_statistics,
                        constructor_statistics=constructor_statistics,
                        execution_statistics=measured_statistics,
                        constructor_setup_statistics=constructor_setup,
                        setup_statistics=setup)
        # All further work is untimed verification; it can raise whole-worker
        # driver memory but cannot change the above RSS/callback snapshots.
        initial = data.unknowns()
        zero = np.zeros(size)
        arrays = {'initial_unknowns': initial, 'zero_residual': data.residual(zero),
                  'zero_jvp': data.jvp(zero, inputs['direction'])}
        arrays.update({'initial_' + key: value for key, value in first_sample.items()})
        data.set_unknowns(inputs['unknowns'])
        arrays['nonzero_residual'] = data.residual(inputs['unknowns'])
        arrays['nonzero_jvp'] = data.jvp(inputs['unknowns'], inputs['direction'])
        arrays.update(data.sample_with_derivatives(inputs['points']))
        work = data.work_statistics()
        diagnostics = data.diagnostics()
        np.savez(args.state_output, **arrays)
    if digest(args.input) != config_hash or digest(args.input_arrays) != arrays_hash:
        raise ValueError('frozen worker inputs changed after decoding')
    result = dict(snapshot, config=values, geometry=args.geometry, execution='kokkos',
                  compiled_execution=actual_space, cpu_threads=actual_threads,
                  requested_threads=args.threads, device=device,
                  library_sha256=backend.library_sha256(),
                  dependency_images=backend.dependency_images,
                  runtime_images=loaded_kokkos_images(), loaded_image_verified=True,
                  input_sha256=config_hash, input_arrays_sha256=arrays_hash,
                  all_arrays_finite=all(np.all(np.isfinite(v)) for v in arrays.values()),
                  initial_unknowns_zero=bool(np.array_equal(initial, zero)),
                  work_statistics=work, diagnostics=diagnostics,
                  memory_window='RSS is lifetime high-water through first sample; Kokkos counters reset after initialization; both captured before verification',
                  coefficient_timing_is_subset_of_first_sample=True,
                  note='No nonlinear solve. Callback peaks exclude CUDA stack/runtime; whole-worker driver observations also include verification.')
    write_json(Path(args.worker_output), result)


def state_shapes(grid):
    shapes = {key: (4*int(np.prod(grid)),) for key in
              ('initial_unknowns', 'zero_residual', 'zero_jvp', 'nonzero_residual', 'nonzero_jvp')}
    for prefix in ('', 'initial_'):
        for key, _ in Point._fields_:
            shape = (len(SAMPLE_POINTS), 9) if key in ('gamma', 'Kij', 'conformal_metric', 'Atilde') else \
                    (len(SAMPLE_POINTS), 4) if key == 'correction' else (len(SAMPLE_POINTS),)
            shapes[prefix+key] = shape
    shapes['dgamma'] = (len(SAMPLE_POINTS), 3, 3, 3)
    return shapes


def compare_arrays(host, execution, grid):
    result = {}
    shapes = state_shapes(grid)
    with np.load(host, allow_pickle=False) as a, np.load(execution, allow_pickle=False) as b:
        if set(a.files) != set(shapes) or set(b.files) != set(shapes):
            return dict(passed=False, reason='declared array keys missing or unexpected', arrays={})
        for key in a.files:
            x, y = a[key], b[key]
            finite = bool(np.all(np.isfinite(x)) and np.all(np.isfinite(y)))
            compatible = x.shape == y.shape == shapes[key]
            dtype = x.dtype == y.dtype == np.dtype('float64')
            difference = float(np.max(np.abs(x-y)/(1+np.abs(x)))) if finite and compatible else None
            result[key] = dict(finite=finite, declared_shape=compatible, dtype_float64=dtype,
                               max_scaled_difference=difference, passed=finite and compatible and dtype and difference <= BOUND)
    return dict(bound=BOUND, arrays=result, passed=bool(result) and all(r['passed'] for r in result.values()))


def protocol_checks(record, variant, images, expected_config, input_hash, arrays_hash):
    setup = record.get('setup_statistics') or {}
    primary = str(Path(variant['hispid_library']).resolve())
    runtime = {str(Path(path).resolve()): images[str(Path(path).resolve())]
               for path in variant.get('runtime_images', [])}
    checks = dict(configuration=record['config'] == expected_config,
                  inputs=record['input_sha256'] == input_hash and record['input_arrays_sha256'] == arrays_hash,
                  execution=record['execution'] == 'kokkos' and record['compiled_execution'] == variant['space']
                            and record['cpu_threads'] == variant['threads'],
                  geometry=setup.get('geometry_execution') == int(record['geometry'] == 'execution'),
                  precision=setup.get('scalar_digits') == 53 if record['geometry'] == 'execution' else setup.get('scalar_digits', 0) >= 53,
                  image=record['loaded_image_verified'] and record['library_sha256'] == images[primary]
                        and all(images.get(p) == sha for p, sha in record['dependency_images'].items()),
                  runtime=record['runtime_images'] == runtime,
                  finite=record['all_arrays_finite'], zero=record['initial_unknowns_zero'],
                  no_solve=record['diagnostics']['newton_iterations'] == 0
                           and record['diagnostics']['krylov_iterations'] == 0,
                  memory_tracking=record['execution_statistics']['memory_tracking_available'] == 1)
    times = [record[k] for k in ('initialization_seconds', 'creation_seconds', 'first_sample_seconds', 'ready_to_sample_seconds')]
    times += [setup.get(k, -1) for k in ('spectral_seconds', 'geometry_seconds', 'coefficient_seconds')]
    checks['timers'] = all(np.isfinite(v) and v >= 0 for v in times)
    if variant['space'] == 'Cuda':
        device = record.get('device') or {}
        checks['one_gpu'] = device.get('visible_count') == 1 and device.get('visible_ordinal') == 0 and bool(device.get('uuid'))
        observed = {r['uuid'] for r in record['driver_memory']['samples']}
        checks['driver_identity'] = not observed or observed == {device.get('uuid')}
    checks['passed'] = all(checks.values())
    return checks


def setup_control(variant, images, raw, timeout, env):
    """Execute both prerequisite suites and bind their executables and DSOs."""
    build = Path(variant['hispid_library']).resolve().parent
    log = raw / (variant['id'] + '_setup_control.log')
    junit = raw / (variant['id'] + '_setup_control.xml')
    leftovers = [p for p in (log, junit) if p.exists()]
    if leftovers:
        orphan = raw / ('interrupted_control_' + str(time.time_ns()) + '_' + variant['id'])
        orphan.mkdir()
        for path in leftovers:
            path.rename(orphan/path.name)
    listing = subprocess.check_output(['ctest', '--test-dir', str(build), '--show-only=json-v1'], text=True)
    tests = {t['name']: t for t in json.loads(listing)['tests'] if t['name'] in SETUP_TESTS}
    if set(tests) != SETUP_TESTS:
        raise ValueError('required native setup/seed suites missing')
    executables = {str(Path(t['command'][0]).resolve(strict=True)): digest(t['command'][0]) for t in tests.values()}
    # Check each test's resolved HiSpID/TwoPunctures/Kokkos dependencies too.
    # This prevents a same-named test from silently loading another build.
    for executable in executables:
        listing = subprocess.check_output(['ldd', executable], text=True)
        if 'not found' in listing:
            raise ValueError('unresolved native test dependencies')
        import re
        dependencies = [str(Path(p).resolve(strict=True)) for p in re.findall(r'(?:=>\s+|^\s*)(/[^\s]+)\s+\(', listing, re.M)]
        required = {p: images[p] for p in images if Path(p).name.startswith(('libHiSpID', 'libTwoPunctures', 'libkokkos'))}
        actual = {p: digest(p) for p in dependencies if Path(p).name.startswith(('libHiSpID', 'libTwoPunctures', 'libkokkos'))}
        if actual != required:
            raise ValueError('native test dependency images differ from workers')
    command = ['ctest', '--test-dir', str(build), '-R', '^(hispid_kokkos_setup|hispid_execution_seed_export)$',
               '--output-on-failure', '--output-junit', str(junit), '-j1']
    started = time.time()
    timed_out = False
    with log.open('w') as stream:
        process = subprocess.Popen(command, env=env, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            code = process.wait()
            timed_out = True
    tested = {}
    if junit.exists():
        for test in ET.parse(junit).getroot().iter('testcase'):
            tested[test.attrib['name']] = all(test.find(key) is None for key in ('failure', 'error', 'skipped'))
    row = dict(executed=True, passed=code == 0 and set(tested) == SETUP_TESTS and all(tested.values()),
               tests=tested, exit_code=code, timeout=timed_out, command=command,
               images=images, executables_sha256=executables,
               log=str(log.relative_to(ROOT)), log_sha256=digest(log), elapsed_seconds=time.time()-started)
    if junit.exists():
        row.update(junit=str(junit.relative_to(ROOT)), junit_sha256=digest(junit))
    return row


def verify_sources(result):
    for path, sha in result['binding']['sources_sha256'].items():
        if digest(ROOT/path) != sha:
            raise ValueError('supplemental source changed: ' + path)


def verify_controls(result):
    for row in result['setup_controls'].values():
        for key in ('log', 'junit'):
            if key in row and digest(ROOT / row[key]) != row[key + '_sha256']:
                raise ValueError('native setup evidence changed')
        for path, sha in row.get('executables_sha256', {}).items():
            if digest(path) != sha:
                raise ValueError('native setup executable changed')
        for path, sha in row.get('images', {}).items():
            if digest(path) != sha:
                raise ValueError('native setup image changed')


def control_passed(control, images):
    return (control.get('executed') is True and control.get('passed') is True
            and control.get('exit_code') == 0 and control.get('timeout') is False
            and set(control.get('tests', {})) == SETUP_TESTS and all(control['tests'].values())
            and len(control.get('executables_sha256', {})) == 2 and control.get('images') == images
            and all(control.get(key) and len(control.get(key+'_sha256', '')) == 64 for key in ('log', 'junit')))


def finalize(result, raw, verify=True):
    result['full_artifact_verification'] = False
    result['declared_setup_completed'] = False
    result['all_setup_checks_passed'] = False
    if verify:
        verify_artifacts(result)
        verify_sources(result)
        verify_controls(result)
        result['full_artifact_verification'] = True
    variants = {v['id']: v for v in result['manifest']['variants']}
    comparisons = {}
    # Both sides must have their current protocol checks before pairing.
    for label, record in result['records'].items():
        grid_id = '_'.join(map(str, record['grid']))
        frozen = result['inputs'][grid_id]
        record['checks'] = protocol_checks(record, variants[record['variant']], result['images'][record['variant']],
                                          frozen['config'], frozen['config_sha256'], frozen['arrays_sha256'])
    for label, record in result['records'].items():
        grid_id = '_'.join(map(str, record['grid']))
        prefix = f'{grid_id}_{record["variant"]}_{record["repeat"]}'
        host = result['records'].get(prefix + '_host')
        execution = result['records'].get(prefix + '_execution')
        if host is None or execution is None or prefix in comparisons:
            continue
        control = result['setup_controls'][record['variant']]
        comparison_binding = dict(host_state=host['state_sha256'], execution_state=execution['state_sha256'],
                                  host_checks=host['checks'], execution_checks=execution['checks'],
                                  control_sha256=hashlib.sha256(json.dumps(control, sort_keys=True).encode()).hexdigest())
        cached = result.get('comparisons', {}).get(prefix, {})
        if cached.get('binding') == comparison_binding:
            comparison = dict(cached)
        else:
            comparison = compare_arrays(ROOT / host['state'], ROOT / execution['state'], record['grid'])
        comparison['binding'] = comparison_binding
        qualified = result['full_artifact_verification'] and comparison['passed'] and host['checks']['passed'] and execution['checks']['passed'] and control_passed(control, result['images'][record['variant']])
        comparison['qualified'] = qualified
        comparison.pop('creation_speedup', None)
        comparison.pop('ready_speedup', None)
        if qualified:
            comparison['creation_speedup'] = host['creation_seconds'] / execution['creation_seconds']
            comparison['ready_speedup'] = host['ready_to_sample_seconds'] / execution['ready_to_sample_seconds']
        comparisons[prefix] = comparison
    result['comparisons'] = comparisons
    result['completed_workers'] = len(result['records']) + len(result['failures'])
    result['coverage_complete'] = set(result['records']) | set(result['failures']) == set(result['expected_worker_ids'])
    result['declared_scope'] = result['binding']['grids'] == GRIDS and result['binding']['repeats'] == 3 and set(variants) == VARIANTS
    result['declared_setup_completed'] = result['full_artifact_verification'] and result['declared_scope'] and result['coverage_complete']
    result['all_setup_checks_passed'] = result['full_artifact_verification'] and result['coverage_complete'] and not result['failures'] and len(comparisons)*2 == len(result['expected_worker_ids']) and all(c['qualified'] for c in comparisons.values())


def coordinator(args):
    manifest = json.loads(Path(args.manifest).read_text())
    manifest['variants'] = [v for v in manifest['variants'] if v['id'] != 'reference']
    variants = manifest['variants']
    if not variants or len({v['id'] for v in variants}) != len(variants) or any(v['execution'] != 'kokkos' for v in variants):
        raise ValueError('unique Kokkos variants required')
    grids = [list(map(int, g.split(':'))) for g in args.grids.split(',')]
    if not grids or len({tuple(g) for g in grids}) != len(grids) or any(len(g) != 3 or min(g) < 4 or max(g) > 256 or g[2] % 2 for g in grids):
        raise ValueError('unique grids4..256 with even phi required')
    if args.repeats < 1 or not np.isfinite(args.timeout) or args.timeout <= 0:
        raise ValueError('positive repeats and timeout required')
    out = Path(args.output).resolve()
    raw = ROOT / 'validation/raw' / out.stem
    if (out.exists() or raw.exists()) and not args.resume:
        raise FileExistsError('existing supplemental campaign; use --resume')
    raw.mkdir(parents=True, exist_ok=True)
    images = {v['id']: manifest_images(v) for v in variants}
    sources = [Path(__file__), ROOT/'validation/benchmark_kokkos.py', ROOT/'validation/benchmark_bowen_york.py', ROOT/'validation/benchmark_krylov_matrix.py', ROOT/'validation/benchmark_common_stopping.py']
    sources += list((ROOT/'python').glob('*.py')) + list((ROOT/'examples').glob('*.py'))
    binding = dict(grids=grids, repeats=args.repeats, timeout=args.timeout, bound=BOUND,
                   sources_sha256={str(p.relative_to(ROOT)): digest(p) for p in sources},
                   geometry_choices=['host', 'execution'], nonlinear_solve=False)
    hw = hardware()
    compatible = hardware_class(hw)
    result = json.loads(out.read_text()) if out.exists() else dict(manifest=manifest, manifest_sha256=digest(args.manifest),
             binding=binding, images=images, hardware_class=compatible, allocation_epochs=[], records={}, failures={},
             inputs={}, input_sha256={}, setup_controls={}, comparisons={})
    if result['binding'] != binding or result['images'] != images or result['manifest_sha256'] != digest(args.manifest) or result['hardware_class'] != compatible:
        raise ValueError('resume sources/images/protocol/hardware class changed')
    result['allocation_epochs'].append(hw)
    epoch = len(result['allocation_epochs']) - 1
    expected = [f'{"_".join(map(str,g))}_{v["id"]}_{repeat}_{geometry}' for g in grids
                for repeat in range(args.repeats) for v in variants for geometry in ('host', 'execution')]
    result['expected_worker_ids'] = expected
    if (set(result['records']) | set(result['failures'])) - set(expected):
        raise ValueError('unexpected retained workers')
    # Only the coordinator loads a single image to obtain default config values;
    # every timed worker loads its own image in a fresh process.
    backend = Backend(variants[0]['hispid_library'])
    for grid in grids:
        key = '_'.join(map(str, grid))
        config = configuration(backend, 'moderate', grid, 1e-12)
        config.memory_limit_mib = 32768
        values = as_dict(config)
        config_path, arrays_path = raw / ('input_' + key + '.json'), raw / ('input_' + key + '.npz')
        if key not in result['inputs']:
            write_json(config_path, values)
            np.savez(arrays_path, **frozen_arrays(grid))
            result['inputs'][key] = dict(config=values, config_path=str(config_path.relative_to(ROOT)),
                                        arrays_path=str(arrays_path.relative_to(ROOT)),
                                        config_sha256=digest(config_path), arrays_sha256=digest(arrays_path))
            result['input_sha256'].update({str(p.relative_to(ROOT)): digest(p) for p in (config_path, arrays_path)})
        elif result['inputs'][key]['config'] != values:
            raise ValueError('resume physical inputs changed')
    del backend
    finalize(result, raw)
    write_json(out, result)
    for grid in grids:
        key = '_'.join(map(str, grid))
        for repeat in range(args.repeats):
            for variant in (variants if repeat % 2 == 0 else list(reversed(variants))):
                env = os.environ.copy()
                env.update(OMP_NUM_THREADS=str(variant['threads']), OMP_PROC_BIND='close', OMP_PLACES='cores', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
                env.update(variant.get('environment', {}))
                if variant['id'] not in result['setup_controls']:
                    verify_sources(result)
                    if manifest_images(variant) != images[variant['id']]:
                        raise ValueError('native control images changed')
                    result['setup_controls'][variant['id']] = setup_control(variant, images[variant['id']], raw, args.timeout, env)
                    write_json(out, result)
                for geometry in (('host', 'execution') if repeat % 2 == 0 else ('execution', 'host')):
                    label = f'{key}_{variant["id"]}_{repeat}_{geometry}'
                    if label in result['records'] or label in result['failures']:
                        continue
                    if (raw / 'STOP_AFTER_WORKER').exists():
                        finalize(result, raw)
                        write_json(out, result)
                        print('Checkpointed', result['completed_workers'], 'of', len(expected), flush=True)
                        return
                    if manifest_images(variant) != images[variant['id']]:
                        raise ValueError('worker images changed')
                    verify_sources(result)
                    paths = {k: raw / (label + suffix) for k, suffix in [('worker', '.json'), ('state', '.npz'), ('log', '.log')]}
                    leftovers = [p for p in paths.values() if p.exists()]
                    if leftovers:
                        orphan = raw / ('interrupted_' + str(time.time_ns()) + '_' + label)
                        orphan.mkdir()
                        for p in leftovers:
                            p.rename(orphan / p.name)
                    frozen = result['inputs'][key]
                    command = [sys.executable, str(Path(__file__).resolve()), '--library', variant['hispid_library'],
                               '--input', str(ROOT / frozen['config_path']), '--input-arrays', str(ROOT / frozen['arrays_path']),
                               '--geometry', geometry, '--threads', str(variant['threads']),
                               '--worker-output', str(paths['worker']), '--state-output', str(paths['state'])]
                    print('Starting', label, flush=True)
                    started = time.time()
                    timed_out = False
                    with paths['log'].open('w') as stream:
                        process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
                        telemetry = GPUProcessMemory(process.pid, variant['space'] == 'Cuda')
                        telemetry.start()
                        try:
                            code = process.wait(timeout=args.timeout)
                        except subprocess.TimeoutExpired:
                            os.killpg(process.pid, signal.SIGKILL)
                            code = process.wait()
                            timed_out = True
                        memory = telemetry.finish()
                    row = dict(command=command, allocation_epoch=epoch, driver_memory=memory, elapsed_seconds=time.time()-started,
                               variant=variant['id'], geometry=geometry, repeat=repeat, grid=grid)
                    for k, p in paths.items():
                        if p.exists():
                            row.update({k: str(p.relative_to(ROOT)), k+'_sha256': digest(p)})
                    if code or not paths['worker'].exists() or not paths['state'].exists():
                        result['failures'][label] = dict(row, exit_code=code, timeout=timed_out)
                    else:
                        result['records'][label] = dict(json.loads(paths['worker'].read_text()), **row)
                    finalize(result, raw, verify=False)
                    write_json(out, result)
                    print('Retained', label, 'attempt', result['completed_workers'], 'of', len(expected), flush=True)
    finalize(result, raw)
    write_json(out, result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest'); p.add_argument('--output'); p.add_argument('--resume', action='store_true')
    p.add_argument('--grids', default=','.join(':'.join(map(str,g)) for g in GRIDS)); p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--timeout', type=float, default=3600)
    p.add_argument('--library'); p.add_argument('--input'); p.add_argument('--input-arrays')
    p.add_argument('--geometry', choices=['host', 'execution']); p.add_argument('--threads', type=int, default=1)
    p.add_argument('--worker-output'); p.add_argument('--state-output')
    args = p.parse_args()
    if not 1 <= args.threads <= 16:
        p.error('threads must be1..16')
    if args.library:
        if not all((args.input, args.input_arrays, args.geometry, args.worker_output, args.state_output)):
            p.error('worker requires frozen inputs, geometry and output paths')
        worker(args)
    else:
        if not args.manifest or not args.output:
            p.error('coordinator requires manifest and output')
        coordinator(args)


if __name__ == '__main__':
    main()
