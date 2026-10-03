"""Read and requalify the separate geometry-setup dataset for the LaTeX report."""
import copy
import hashlib
import json
from pathlib import Path
import statistics

import numpy as np

from benchmark_bowen_york import digest
from benchmark_geometry_setup import BOUND, GRIDS, VARIANTS, decode_config, finalize, frozen_arrays
from benchmark_kokkos import manifest_images
from kokkos_build_manifest import native_source, setup_build_source

SPACES = {'serial': ('Serial', 1), 'cuda': ('Cuda', 16),
          **{f'openmp{n}': ('OpenMP', n) for n in (1, 2, 4, 8, 16)}}


def inventory(result):
    """Reconstruct coverage rather than trusting cached completion flags."""
    binding = result['binding']
    if (binding['grids'] != GRIDS or binding['repeats'] != 3
            or binding['geometry_choices'] != ['host', 'execution']
            or binding['nonlinear_solve'] is not False or binding['bound'] != BOUND
            or not np.isfinite(binding['timeout']) or binding['timeout'] <= 0):
        raise ValueError('setup report requires the declared paired 126-worker protocol')
    variants = result['manifest']['variants']
    if len(variants) != len(VARIANTS) or {v['id'] for v in variants} != VARIANTS:
        raise ValueError('setup variants differ from the declared protocol')
    for variant in variants:
        if ((variant['space'], variant['threads']) != SPACES[variant['id']]
                or variant['execution'] != 'kokkos' or not variant.get('runtime_images')):
            raise ValueError('setup execution space/thread count is inconsistent')
    expected = {f'{"_".join(map(str,g))}_{v["id"]}_{repeat}_{geometry}'
                for g in GRIDS for v in variants for repeat in range(3)
                for geometry in ('host', 'execution')}
    declared = result['expected_worker_ids']
    if len(declared) != len(expected) or set(declared) != expected:
        raise ValueError('setup worker inventory differs from the declared protocol')
    records, failures = result['records'], result['failures']
    if set(records) & set(failures) or (set(records) | set(failures)) - expected:
        raise ValueError('setup worker coverage is overlapping or unexpected')
    for label, row in {**records, **failures}.items():
        identity = f'{"_".join(map(str,row["grid"]))}_{row["variant"]}_{row["repeat"]}_{row["geometry"]}'
        if identity != label:
            raise ValueError('setup worker identity differs from its retained row: ' + label)
        for key in (('log', 'worker', 'state') if label in records else ('log',)):
            if not row.get(key) or len(row.get(key+'_sha256', '')) != 64:
                raise ValueError('setup worker evidence is missing: ' + label)
        if row['variant'] not in result['setup_controls']:
            raise ValueError('setup worker native control is missing: ' + label)
    if set(result['setup_controls']) - VARIANTS:
        raise ValueError('unexpected setup native controls')
    return expected


def verify_manifest(result, root):
    manifest = result['manifest']
    encoded = (json.dumps(manifest, indent=2) + '\n').encode()
    if hashlib.sha256(encoded).hexdigest() != result['manifest_sha256']:
        raise ValueError('setup manifest differs from its retained hash')
    if (manifest.get('scope') != 'supplemental paired geometry setup'
            or manifest.get('solve_campaign_replaced') is not False):
        raise ValueError('setup manifest must remain separate from the solve campaign')
    if set(manifest['builds']) != {'serial', 'openmp', 'cuda'}:
        raise ValueError('setup build inventory must contain serial/openmp/cuda')
    for path, sha in manifest['source_sha256'].items():
        if digest(root/path) != sha:
            raise ValueError('setup manifest source changed: ' + path)
    expected_native = native_source(root)
    for build in manifest['builds'].values():
        cache = Path(build['cache_path'])
        if digest(cache) != build['cache_sha256'] or cache.read_text() != build['cache']:
            raise ValueError('setup build cache changed')
        if setup_build_source(build, expected_native) != build['compiled_source']:
            raise ValueError('setup compiled source witness changed')
    if set(result['images']) != VARIANTS:
        raise ValueError('setup image inventory differs from the variants')
    for variant in manifest['variants']:
        if manifest_images(variant) != result['images'][variant['id']]:
            raise ValueError('setup native/runtime images changed: ' + variant['id'])
        build = manifest['builds']['openmp' if variant['id'].startswith('openmp') else variant['id']]
        native = {str(Path(variant[key]).resolve()): result['images'][variant['id']][str(Path(variant[key]).resolve())]
                  for key in ('hispid_library', 'by_library')}
        runtime = {str(Path(path).resolve()): result['images'][variant['id']][str(Path(path).resolve())]
                   for path in variant['runtime_images']}
        if native != build['images'] or runtime != build['runtime_images']:
            raise ValueError('setup variant image bank differs from its compiled build')


def verify_inputs_and_workers(result, root):
    expected_inputs = {'_'.join(map(str,g)) for g in GRIDS}
    if set(result['inputs']) != expected_inputs:
        raise ValueError('setup frozen input inventory is incomplete')
    witnesses = {}
    for grid in GRIDS:
        saved = result['inputs']['_'.join(map(str, grid))]
        for key in ('config', 'arrays'):
            path, sha = saved[key+'_path'], saved[key+'_sha256']
            if digest(root/path) != sha:
                raise ValueError('setup frozen input changed: ' + path)
            witnesses[path] = sha
        values = json.loads((root/saved['config_path']).read_bytes())
        if values != saved['config'] or values['n'] != grid:
            raise ValueError('setup frozen configuration differs from the dataset')
        decode_config(values)
        expected = frozen_arrays(grid)
        with np.load(root/saved['arrays_path'], allow_pickle=False) as arrays:
            if set(arrays.files) != set(expected) or any(
                    arrays[key].dtype != np.dtype('float64')
                    or not np.array_equal(arrays[key], value) for key, value in expected.items()):
                raise ValueError('setup frozen arrays differ from the declared recipe')
    if witnesses != result['input_sha256']:
        raise ValueError('setup input hashes differ from their artifact inventory')
    for label, row in result['records'].items():
        if digest(root/row['worker']) != row['worker_sha256']:
            raise ValueError('setup worker JSON changed: ' + label)
        worker = json.loads((root/row['worker']).read_bytes())
        if not worker or any(key not in row or row[key] != value for key, value in worker.items()):
            raise ValueError('setup record differs from its original worker JSON: ' + label)


def load_setup_results(path, artifact_root):
    """Read-only validation; recompute all 22 array checks and qualifications."""
    path, root = Path(path), Path(artifact_root).resolve(strict=True)
    payload = path.read_bytes()
    sha = hashlib.sha256(payload).hexdigest()
    result = copy.deepcopy(json.loads(payload))
    inventory(result)
    verify_manifest(result, root)
    verify_inputs_and_workers(result, root)
    result['comparisons'] = {}  # Never admit a cached/forged array-comparison flag.
    finalize(result, root, artifact_root=root)
    if digest(path) != sha:
        raise RuntimeError('setup JSON changed during report generation; use an immutable snapshot')
    return result, sha


def qualified_group(result, grid, variant):
    prefix = f'{"_".join(map(str,grid))}_{variant}_'
    pairs = [result['comparisons'].get(prefix+str(repeat), {}) for repeat in range(3)]
    if not result.get('full_artifact_verification') or not all(p.get('qualified') for p in pairs):
        return None
    host = [result['records'][prefix+str(repeat)+'_host'] for repeat in range(3)]
    execution = [result['records'][prefix+str(repeat)+'_execution'] for repeat in range(3)]
    for key in ('creation_seconds', 'ready_to_sample_seconds'):
        if any(row[key] <= 0 for row in host+execution):
            return None
    med = lambda rows, key: statistics.median(row[key] for row in rows)
    return dict(creation_speedup=med(host, 'creation_seconds')/med(execution, 'creation_seconds'),
                ready_speedup=med(host, 'ready_to_sample_seconds')/med(execution, 'ready_to_sample_seconds'),
                rss_saved_percent=100*(1-med(execution, 'max_rss_bytes')/med(host, 'max_rss_bytes')))


def render_setup(result, sha, tex, number):
    doc = [r'''\section{Separate execution-space geometry setup}
This supplemental campaign uses newer images and varies only host versus
execution-space geometry construction within each image. It measures spinning
and boosted seed geometry, coordinates, differentiation and coefficient caches,
construction and first sampling. It performs no nonlinear solve. Its results do
not replace the frozen solve matrix or qualify an extreme physical binary.
The declared matrix has 126 workers: three grids, seven execution/thread
variants, two geometry choices and three matched repeats. Both native setup
and public seed export suites must execute without a skipped wider-precision
oracle. All 22 retained arrays are reread, with exact declared shapes,
float64 storage, finite values and a $10^{-10}$ scaled-difference bound.
Cached completion and comparison flags are not used to establish qualification.

RSS is the lifetime high-water through first sample, including input loading.
Kokkos callback peaks reset after initialization and retain any tracked
allocations still live at reset; native vectors and CUDA stack/runtime outside
the callbacks are excluded. Driver observations cover the entire worker,
including later verification, at 1 Hz, and are lower bounds on transient peaks.
Ready time adds construction and first sampling; coefficient transformation is
a subset of first sampling. The metadata-query gap is reported separately.
Physical sampling still includes CPU work. Native block factorization belongs
to the separate nonlinear solve campaign and is absent from these setup timings.
''']
    attempts = len(result['records'])+len(result['failures'])
    bad = sum(not row['checks']['passed'] for row in result['records'].values())
    rejected = sum(not pair['qualified'] for pair in result['comparisons'].values())
    doc.append(f'Retained setup attempts: {attempts}/126; process failures: {len(result["failures"])}; '
               f'failed worker checks: {bad}; unqualified retained pairs: {rejected}.'+r'\par'+'\n')
    if not result['declared_setup_completed']:
        doc.append(r'\textbf{Setup campaign pending. Missing workers are not estimated.}\par'+'\n')
    doc.append(r'\small\begin{longtable}{llrrrr}\toprule Grid & Variant & Pairs & Create speedup & Ready speedup & RSS saved\\\midrule\endhead'+'\n')
    for grid in GRIDS:
        for variant in result['manifest']['variants']:
            label = variant['id']
            pairs = sum(f'{"_".join(map(str,grid))}_{label}_{repeat}' in result['comparisons'] for repeat in range(3))
            ratio = qualified_group(result, grid, label)
            cells = [r'$\times$'.join(map(str,grid)), tex(label), str(pairs),
                     number(ratio['creation_speedup']) if ratio else '--',
                     number(ratio['ready_speedup']) if ratio else '--',
                     number(ratio['rss_saved_percent'])+r'\%' if ratio else '--']
            doc.append(' & '.join(cells)+r'\\'+'\n')
    doc.append(r'\bottomrule\end{longtable}\normalsize'+'\n')
    doc.append('Ratios require all three matched pairs to qualify. Negative RSS savings mean increased host memory. '
               'All other timings are diagnostic; a dagger marks a group with a failed worker check.\n')
    for grid in GRIDS:
        rows = [row for row in result['records'].values() if row['grid'] == grid]
        if not rows:
            continue
        doc.append(r'\subsection{Setup grid '+r'$\times$'.join(map(str,grid))+'}\n')
        doc.append(r'\scriptsize\begin{longtable}{llrrrrrrr}\toprule Variant & Geometry & $n$ & Init(s) & Create(s) & Sample(s) & Ready(s) & RSS(GiB) & Driver(GiB)\\\midrule\endhead'+'\n')
        groups = []
        for variant in result['manifest']['variants']:
            for geometry in ('host', 'execution'):
                group = [r for r in rows if r['variant'] == variant['id'] and r['geometry'] == geometry]
                if not group:
                    continue
                groups.append((variant['id'], geometry, group))
                med = lambda key: statistics.median(row[key] for row in group)
                driver = [r['driver_memory']['sampled_process_peak_bytes']/2**30 for r in group
                          if r['driver_memory'].get('sampled_process_peak_bytes') is not None]
                label = tex(variant['id'])+(r'$\dagger$' if not all(r['checks']['passed'] for r in group) else '')
                cells = [label, geometry, str(len(group)), *[number(med(k)) for k in
                         ('initialization_seconds', 'creation_seconds', 'first_sample_seconds', 'ready_to_sample_seconds')],
                         number(med('max_rss_bytes')/2**30), number(statistics.median(driver)) if driver else '--']
                doc.append(' & '.join(cells)+r'\\'+'\n')
        doc.append(r'\bottomrule\end{longtable}\normalsize'+'\n')
        doc.append(r'\scriptsize\begin{longtable}{llrrrrrrr}\toprule Variant & Geometry & Spectral(s) & Geometry(s) & Coeff(s) & Gap(s) & Wall(s) & K host(GiB) & K device(GiB)\\\midrule\endhead'+'\n')
        for variant, geometry, group in groups:
            phases = [number(statistics.median(row['setup_statistics'][key] for row in group))
                      for key in ('spectral_seconds', 'geometry_seconds', 'coefficient_seconds')]
            wall = [number(statistics.median(row[key] for row in group))
                    for key in ('instrumentation_gap_seconds', 'wall_creation_through_sample_seconds')]
            memory = [number(statistics.median(row['execution_statistics'][key]/2**30 for row in group))
                      for key in ('kokkos_host_peak_bytes', 'kokkos_device_peak_bytes')]
            doc.append(' & '.join([tex(variant), geometry, *phases, *wall, *memory])+r'\\'+'\n')
        doc.append(r'\bottomrule\end{longtable}\normalsize'+'\n')
    for variant, control in result['setup_controls'].items():
        doc.append('Native suites '+tex(variant)+': '+('passed' if control['passed'] else 'failed or skipped')+r'.\par'+'\n')
    if result['failures']:
        doc.append(r'\small\begin{longtable}{p{0.6\linewidth}lrr}\toprule Setup worker & Outcome & Elapsed(s) & Exit\\\midrule\endhead'+'\n')
        for label, row in result['failures'].items():
            doc.append(' & '.join([r'\nolinkurl{'+label+'}', 'timeout' if row.get('timeout') else 'failed',
                       number(row.get('elapsed_seconds')), tex(row.get('exit_code', '--'))])+r'\\'+'\n')
        doc.append(r'\bottomrule\end{longtable}\normalsize'+'\n')
    doc.append('Setup JSON SHA256:'+r'\par{\footnotesize\ttfamily '+sha+r'}\par'+'\n')
    doc.append('Setup manifest SHA256:'+r'\par{\footnotesize\ttfamily '+result['manifest_sha256']+r'}\par'+'\n')
    for variant in result['manifest']['variants']:
        for path in (variant['hispid_library'], variant['by_library'], *variant['runtime_images']):
            key = str(Path(path).resolve())
            doc.append(tex(variant['id'])+': '+r'\nolinkurl{'+key+r'}\par{\footnotesize\ttfamily '+result['images'][variant['id']][key]+r'}\par'+'\n')
    return ''.join(doc)
