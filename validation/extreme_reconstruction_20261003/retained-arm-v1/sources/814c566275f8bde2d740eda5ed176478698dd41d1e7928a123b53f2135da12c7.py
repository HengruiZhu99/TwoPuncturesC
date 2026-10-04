"""Bound, fixed-state Cartesian reconstruction diagnosis; never re-solves data."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np
from checkpoint_export import read_checkpoint
from configs import as_dict


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('checkpoint', 'solve-artifact', 'point-field-artifact', 'parent-result',
                'binding-audit', 'row-receipt', 'executable', 'build-receipt', 'output'):
        p.add_argument('--' + key, required=True)
    a = p.parse_args()
    paths = {key: Path(getattr(a, key.replace('-', '_'))).resolve(strict=True)
             for key in ('checkpoint', 'solve-artifact', 'point-field-artifact', 'parent-result',
                         'binding-audit', 'row-receipt', 'executable', 'build-receipt')}
    parent = json.loads(paths['parent-result'].read_text())
    audit = json.loads(paths['binding-audit'].read_text())
    row = json.loads(paths['row-receipt'].read_text())['record']
    build = json.loads(paths['build-receipt'].read_text())
    cfg, unknowns, cp = read_checkpoint(paths['checkpoint'])
    if (not parent['completed'] or parent['solve_performed'] or audit.get('passed') is not True
            or audit['parent_result_sha256'] != digest(paths['parent-result'])
            or audit['checkpoint_sha256'] != cp['file_sha256']
            or audit['row_receipt_sha256'] != digest(paths['row-receipt'])
            or cp['source_library_sha256'] != row['library_sha256']
            or cp['parameterization'] != row['unknown_parameterization_id']
            or cp['acceptance'] != 'diagnostic' or as_dict(cfg) != row['config']):
        raise ValueError('failed iterate, parent replay and binding audit must agree')
    expected_solve = next(sha for path, sha in row['raw_artifact_sha256'].items() if path.endswith('_solve.npz'))
    if digest(paths['solve-artifact']) != expected_solve:
        raise ValueError('retained solve artifact changed')
    with np.load(paths['solve-artifact']) as raw:
        saved = raw['unknowns']
        if unknowns.shape != saved.shape or unknowns.dtype != saved.dtype or not np.array_equal(
                unknowns.view(np.uint64), saved.view(np.uint64)):
            raise ValueError('retained coefficient bits differ')
    fd = next(r for r in parent['FD_records'] if r['step_factor'] == 1)
    if digest(paths['point-field-artifact']) != fd['raw_sha256']:
        raise ValueError('original sampled point fields changed')
    if digest(paths['executable']) != build['executable_sha256']:
        raise ValueError('diagnostic executable changed')
    bound = {str(path): digest(path) for path in paths.values()}
    bound.update(build['bound_artifacts_sha256'])
    bound[str(Path(__file__).resolve())] = digest(__file__)
    for module in ('checkpoint_export', 'configs', 'hispid'):
        path = Path(__import__(module).__file__).resolve()
        bound[str(path)] = digest(path)

    def verify():
        if any(digest(path) != sha for path, sha in bound.items()):
            raise ValueError('bound diagnosis input, source, image or output changed')

    verify()
    root = Path(a.output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    points = np.asarray(parent['node_xyz'] + parent['worst_offgrid_xyz'])
    nodes = parent['node_indices'] + [-1] * len(parent['worst_offgrid_xyz'])
    with np.load(paths['point-field-artifact']) as raw:
        if not np.array_equal(raw['xyz'], points):
            raise ValueError('parent point inventory differs')
        original_fields = {k: raw[k].copy() for k in ('sampled_gamma', 'sampled_Kij', 'sampled_psi', 'sampled_dgamma', 'H', 'M')}
    point_path = root / 'points.txt'
    point_path.write_text(str(len(points)) + '\n' + ''.join(
        str(node) + ' ' + ' '.join(format(float(x), '.17g') for x in xyz) + '\n'
        for node, xyz in zip(nodes, points)))
    bound[str(point_path)] = digest(point_path)
    output = root / 'Cartesian.json'
    cmd = [str(paths['executable']), str(paths['checkpoint']), cp['source_library_sha256'],
           str(point_path), str(output)]
    result = dict(purpose='fixed_iterate_Cartesian_reconstruction_diagnostic',
                  solve_performed=False, physical_acceptance=False, completed=False,
                  checkpoint=cp, bound_artifacts_sha256=bound, command=cmd,
                  criteria=dict(six_point_sampler_scaled_linf=1e-12),
                  original_result_modified=False)

    def save():
        (root / 'result.json').write_text(json.dumps(result, indent=2) + '\n')

    save()
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    log_path = root / 'run.log'
    with log_path.open('w') as log:
        try:
            run = subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=600)
            result['returncode'] = run.returncode
        except subprocess.TimeoutExpired:
            result['returncode'] = 'timeout'
    bound[str(log_path)] = digest(log_path)
    if output.exists():
        bound[str(output)] = digest(output)
    save()
    verify()
    if result['returncode'] != 0:
        return 1
    numerical = json.loads(output.read_text())
    records = numerical['records']
    if (len(records) != len(points) or numerical['parameterization'] != cp['parameterization']
            or numerical['residual_scaling'] != row['residual_scaling']
            or numerical['solve_performed'] is not False or numerical['physical_acceptance'] is not False
            or not np.array_equal(np.array([r['lab_xyz'] for r in records]), points)):
        raise ValueError('Cartesian diagnostic inventory differs')
    differences = {}
    for key in ('gamma', 'Kij', 'psi', 'dgamma'):
        old = original_fields['sampled_' + key]
        new = np.asarray([r['sampled_' + key] for r in records]).reshape(old.shape)
        if not np.isfinite(new).all():
            raise ValueError('nonfinite reconstructed field')
        differences[key] = float(np.max(np.abs(old - new) / (1 + np.abs(old))))
    result['six_point_sampler_differences'] = differences
    result['six_point_sampler_passed'] = max(differences.values()) <= 1e-12
    result['scalar_digits'] = numerical['scalar_digits']
    result['continuous_tensor_HM'] = [r['continuous']['tensor_constraints'] for r in records]
    result['continuous_equation_HM'] = [r['continuous']['equation_physical'] for r in records]
    result['original_FD_H'] = original_fields['H'].tolist()
    result['original_FD_M'] = original_fields['M'].tolist()
    result['node_derivative_scaled_errors'] = {}
    for label, begin, end in (('value', 0, 1), ('gradient', 1, 4), ('Hessian', 4, 10)):
        old = np.asarray([r['collocation_jets'] for r in records if r['source_node_index'] >= 0])[:, :, begin:end]
        new = np.asarray([r['Cartesian_jets'] for r in records if r['source_node_index'] >= 0])[:, :, begin:end]
        result['node_derivative_scaled_errors'][label] = float(np.max(np.abs(old - new) / (1 + np.abs(old))))
    result['completed'] = True
    result['note'] = ('Diagnostic Cartesian polynomial Hessians bypass collocation differentiation and A/B transforms. '
                      'Native seed/tensor algebra is reused; this is not an independent physical acceptance gate. '
                      'The sampler comparison covers only six retained points and no acceptance transfers.')
    verify()
    save()
    print(json.dumps({key: result[key] for key in ('scalar_digits', 'six_point_sampler_differences',
                       'node_derivative_scaled_errors', 'continuous_tensor_HM', 'original_FD_H')}), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
