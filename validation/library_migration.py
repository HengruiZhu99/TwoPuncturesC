"""Compare byte-bound replay in separate native processes.

dyld can coalesce copied libraries with the same embedded install name. Each
worker verifies its loaded image. The coordinator loads no native library.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
from hispid import Backend
from checkpoints import ROOT, library_sha, restore_payload, select_record
from physical import constraints, norms
from run_validation import points


def worker(library, source_sha, cases, output, axis, equation_case, equation_resolution, equation_nphi):
    def guarded_payload(backend,record):
        token=record.get('unknown_parameterization_id','W_plus_Aminus1_V')
        if token!=backend.parameterization() or record.get('collocation_maps')!=backend.parameterization_maps():
            raise ValueError('different continuous basis/maps cannot be an API migration')
        return restore_payload(record,backend.config())
    backend = Backend(library)
    arrays, reports = {}, []
    for label in cases:
        selector=label.split('@')
        if len(selector) not in (1,3):raise ValueError('case selector is name or name@resolution@nphi')
        record = select_record(selector[0],*(map(int,selector[1:]) if len(selector)==3 else ()))
        if record['library_sha256'] != source_sha:
            raise ValueError('source checkpoint/build mismatch')
        cfg, unknowns = guarded_payload(backend,record)
        xyz, near, bulk = points(cfg, record['horizon_scaled'])
        xyz = xyz[:near+bulk]
        create = backend.create_sampler if hasattr(backend.lib, 'HiSpID_create_sampler') else backend.create
        with create(cfg) as solution:
            solution.set_unknowns(unknowns)
            for name, values in solution.sample(xyz).items(): arrays[label+'__'+name] = values
            details = dict(case=record['case'],array_key=label,resolution=record['resolution'],
                           original_evidence_sha=record['library_sha256'])
            if axis:
                center = .5*(np.array(cfg.hole[0].center)+np.array(cfg.hole[1].center))
                half = .5*(np.array(cfg.hole[0].center)-np.array(cfg.hole[1].center))
                locations = center+np.array([-2., -.7, 0, .7, 2.])[:, None]*half
                checks = []
                for step in (.008, .004, .002, .001):
                    residual = constraints(solution.sample, locations, step)
                    checks.append(dict(step=step, norms=norms(residual),
                                       H=residual['H'].tolist(), M=residual['M'].tolist(),
                                       g=residual['attenuation'].tolist()))
                details.update(axis_points=locations.tolist(), axis_checks=checks)
        reports.append(details)
    record = select_record(equation_case,equation_resolution,equation_nphi)
    if record['library_sha256'] != source_sha: raise ValueError('equation checkpoint/build mismatch')
    cfg, unknowns = guarded_payload(backend,record)
    direction = np.random.default_rng(14108607).normal(0, 1e-4, len(unknowns))
    with backend.create(cfg) as solution:
        arrays['residual'] = solution.residual(unknowns)
        arrays['jvp'] = solution.jvp(unknowns, direction)
    np.savez_compressed(output, **arrays)
    output.with_suffix('.json').write_text(json.dumps(dict(
        loaded_library=str(backend.path), library_sha256=library_sha(backend),
        cases=reports), indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--old-library'); parser.add_argument('--new-library')
    parser.add_argument('--cases', default='moderate_far0_stable,highspin_stable')
    parser.add_argument('--output', default='validation/axis_api_migration.json')
    parser.add_argument('--worker-library'); parser.add_argument('--source-sha')
    parser.add_argument('--equation-case',default='moderate_far0_stable');parser.add_argument('--equation-resolution',type=int,default=24);parser.add_argument('--equation-nphi',type=int,default=12)
    parser.add_argument('--change',default='axis value convention, sampling-only context and physical metric-gradient API')
    parser.add_argument('--equation-atol',type=float,default=0,help='predeclared absolute bound for intentionally changed floating-point evaluation; fields still require bitwise agreement')
    parser.add_argument('--worker-output'); parser.add_argument('--axis', action='store_true')
    args = parser.parse_args(); cases = args.cases.split(',')
    if not np.isfinite(args.equation_atol) or args.equation_atol<0:raise ValueError('nonnegative finite equation bound required')
    if args.worker_library:
        worker(args.worker_library, args.source_sha, cases, Path(args.worker_output), args.axis,args.equation_case,args.equation_resolution,args.equation_nphi)
        return 0
    if not args.old_library or not args.new_library: parser.error('both native libraries are required')
    paths = [Path(args.old_library).resolve(strict=True), Path(args.new_library).resolve(strict=True)]
    hashes = [hashlib.sha256(path.read_bytes()).hexdigest() for path in paths]
    if hashes[0] == hashes[1]: raise ValueError('migration requires distinct binaries')
    temporary = ROOT/'validation/raw'/Path(args.output).stem; temporary.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy(); env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    files = []
    for side, path in zip(('old', 'new'), paths):
        output = temporary/(side+'.npz'); files.append(output)
        command = [sys.executable, str(Path(__file__).resolve()), '--worker-library', str(path),
                   '--source-sha', hashes[0], '--worker-output', str(output), '--cases', args.cases,'--equation-case',args.equation_case,'--equation-resolution',str(args.equation_resolution),'--equation-nphi',str(args.equation_nphi)]
        if side == 'new': command.append('--axis')
        subprocess.run(command, cwd=ROOT, env=env, check=True, timeout=1200)
    reports = [json.loads(path.with_suffix('.json').read_text()) for path in files]
    if [report['library_sha256'] for report in reports] != hashes: raise ValueError('worker build mismatch')
    result = dict(old_library_sha256=hashes[0], new_library_sha256=hashes[1],
                  isolated_processes=True, loaded_images_verified=True,
                  change=args.change,
                  cases=reports[1]['cases'], passed_off_axis_and_equations=False)
    with np.load(files[0]) as old, np.load(files[1]) as new:
        if set(old.files) != set(new.files): raise ValueError('worker array inventory differs')
        for case in result['cases']:
            label = case.get('array_key',case['case'])
            case['off_axis_absolute_differences'] = {
                name.split('__', 1)[1]: float(np.max(abs(new[name]-old[name])))
                for name in old.files if name.startswith(label+'__')}
        result['equation_residual_absolute_difference'] = float(np.max(abs(new['residual']-old['residual'])))
        result['jvp_absolute_difference'] = float(np.max(abs(new['jvp']-old['jvp'])))
    result['equation_comparison_atol']=args.equation_atol
    result['equations_bitwise_identical']=result['equation_residual_absolute_difference']==0 and result['jvp_absolute_difference']==0
    result['passed_off_axis_and_equations'] = (
        all(max(case['off_axis_absolute_differences'].values()) == 0 for case in result['cases'])
        and result['equation_residual_absolute_difference'] <= args.equation_atol and result['jvp_absolute_difference'] <= args.equation_atol)
    (ROOT/args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)
    return 0 if result['passed_off_axis_and_equations'] else 1


if __name__ == '__main__': raise SystemExit(main())
