"""One retained grid of a separate aligned-spin or head-on boost investigation.

Requires the completed performance matrix and a receipt for its successfully
compiled standalone report. Uses the common physical validation machinery with
explicit Kokkos execution. Every export remains diagnostic until the separate
charge/covariance/horizon gates qualify the complete evidence bundle.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from checkpoint_export import write_checkpoint
from configs import as_dict
from execution import select,name,concurrency,device_description
from hispid import Backend,Hole
from run_validation import solve_case
from native_loader import loaded_kokkos_images,validate_krylov


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frozen_inputs(paths):
    payloads={};hashes={}
    for key,path in paths.items():
        data=Path(path).read_bytes();hashes[key]=hashlib.sha256(data).hexdigest()
        if key!='report':payloads[key]=json.loads(data)
    verify_hashes(paths,hashes)
    return payloads,hashes


def verify_hashes(paths,hashes):
    if any(digest(path)!=hashes[key] for key,path in paths.items()):
        raise ValueError('bound prerequisite or attempt evidence changed')


def prerequisites(performance,receipt,hashes):
    if (performance.get('declared_performance_completed') is not True
        or performance.get('expected_workers')!=288 or performance.get('completed_workers')!=288
        or len(performance.get('records',{}))+len(performance.get('failures',{}))!=288):
        raise ValueError('complete declared 288-worker performance matrix required')
    if (receipt.get('compilation_confirmed') is not True
        or receipt.get('compiler')!='mcp__codex_app__compile_latex_document'
        or receipt.get('performance_sha256')!=hashes['performance']
        or receipt.get('report_sha256')!=hashes['report']):
        raise ValueError('matching successful native LaTeX compilation receipt required')
    return performance,receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library',required=True);p.add_argument('--performance-results',required=True)
    p.add_argument('--compiled-report-receipt',required=True);p.add_argument('--compiled-report',required=True)
    p.add_argument('--seed-controls',required=True);p.add_argument('--plan',default=str(Path(__file__).with_name('extreme_kokkos_plan.json')))
    p.add_argument('--case',choices=('aligned_spin99_kokkos','headon_gamma10_kokkos'),required=True)
    p.add_argument('--grid-index',type=int,required=True)
    p.add_argument('--separation',type=float,help='retain a separately labeled separation-calibration case')
    p.add_argument('--output-directory',required=True);p.add_argument('--threads',type=int,default=16)
    p.add_argument('--allow-diagnostic-investigation',action='store_true',help='retain failed prerequisite gates while measuring unqualified data')
    a=p.parse_args()
    input_paths=dict(performance=a.performance_results,receipt=a.compiled_report_receipt,
        report=a.compiled_report,plan=a.plan,seed=a.seed_controls)
    inputs,input_hashes=frozen_inputs(input_paths)
    performance,receipt=prerequisites(inputs['performance'],inputs['receipt'],input_hashes)
    if not 1<=a.threads<=16:raise ValueError('one allocated GPU with at most16 host threads required')
    plan=inputs['plan'];case=next(c for c in plan['cases'] if c['label']==a.case)
    controls=plan['solve_controls'];validate_krylov(controls['krylov'],controls['linear_rtol'])
    grids=case['grids']+([case['independent_fourier_control']] if 'independent_fourier_control' in case else [])
    if not 0<=a.grid_index<len(grids):raise ValueError('grid index outside declared case')
    if a.separation is not None and (not np.isfinite(a.separation) or a.separation<=0):raise ValueError('positive finite separation required')
    library=Path(a.library).resolve(strict=True)
    variant=next(v for v in performance['manifest']['variants'] if v['id']=='cuda')
    measured=performance['images']['cuda']
    if any(digest(path)!=sha for path,sha in measured.items()):raise ValueError('measured CUDA build/dependencies changed')
    if digest(library)!=measured[str(Path(variant['hispid_library']).resolve())]:
        raise ValueError('producer differs from the measured CUDA build')
    seed=inputs['seed']
    if seed.get('library_sha256')!=digest(library) or {c['case'] for c in seed.get('cases',[])}!={'spin99','gamma10'}:
        raise ValueError('fresh bound chi=.99/Gamma=10 seed controls required')
    inputs={c['case']:c for c in seed['cases']}
    if inputs['spin99']['seed_rest_chi']!=.99 or abs(inputs['gamma10']['input_lorentz_factor']/10-1)>1e-13:
        raise ValueError('seed controls do not represent the requested targets')
    prerequisites_passed=bool(seed.get('passed') and performance.get('all_stopping_checks_passed')
                              and performance.get('strict_state_comparisons_passed'))
    if not prerequisites_passed and not a.allow_diagnostic_investigation:
        raise ValueError('failed prerequisites require an explicitly diagnostic investigation')
    os.environ.update(OMP_NUM_THREADS=str(a.threads),OPENBLAS_NUM_THREADS='1',OMP_PROC_BIND='close',OMP_PLACES='cores')
    backend=Backend(str(library));select(backend.lib,'kokkos',a.threads)
    if (backend.parameterization_maps()!=dict(radial_stretch=plan['common_free_data']['maps'][0],angular_stretch=plan['common_free_data']['maps'][1])
        or backend.residual_scaling()!='sin3_alpha_beta' or plan['solve_controls']['row_power']!=3):
        raise ValueError('producer basis/maps/row scaling differs from the declared plan')
    actual_images=backend.dependency_images | loaded_kokkos_images()
    if any(measured.get(path)!=sha for path,sha in actual_images.items()):
        raise ValueError('loaded puncture/runtime image differs from the measured CUDA build')
    device=device_description(backend.lib)
    if name(backend.lib)!='Cuda' or concurrency(backend.lib)!=a.threads or not device or device['visible_count']!=1 or device['visible_ordinal']!=0:
        raise ValueError('actual one-GPU CUDA execution and requested host concurrency required')
    initial_separation=12. if a.case.startswith('aligned') else 25.
    separation=a.separation if a.separation is not None else initial_separation
    label=a.case+(('_d'+format(separation,'.17g').replace('.','p')) if a.separation is not None else '')
    fourier=a.grid_index==len(case['grids'])
    if fourier:label+='_fourier_control'
    root=Path(a.output_directory).resolve();root.mkdir(parents=True,exist_ok=True)
    report_path=root/'results.json';raw=root/'raw';raw.mkdir(exist_ok=True)
    report=json.loads(report_path.read_text()) if report_path.exists() else {}
    previous=report.get(label,{}).get('records',[])
    grid=grids[a.grid_index]
    if any(r['resolution']==grid for r in previous):raise ValueError('grid already retained; use a separate case for a new attempt')
    if any((raw/f'{label}_{grid[0]}_{grid[2]}{suffix}.npz').exists() for suffix in ('','_collocation')):
        raise FileExistsError('unrecorded attempt artifacts exist; preserve them and use a fresh output location')
    if not fourier and a.grid_index!=len(previous):raise ValueError('retain the declared radial sequence in order')
    if [r['resolution'] for r in previous]!=case['grids'][:len(previous)]:
        raise ValueError('retained resolutions differ from the declared grid prefix')
    binding=dict(schema='hispid_extreme_attempt_binding_v1',case=a.case,label=label,
        prerequisite_sha256=input_hashes,measured_images=measured,loaded_images=actual_images,
        producer=str(library),coordinate_separation=separation,host_threads=a.threads)
    binding_path=root/(label+'_attempt_binding.json')
    if binding_path.exists():
        bound_bytes=binding_path.read_bytes()
        if json.loads(bound_bytes)!=binding or digest(binding_path)!=hashlib.sha256(bound_bytes).hexdigest():
            raise ValueError('retained attempt binding differs; use a fresh output location')
    else:
        if previous:raise ValueError('retained rows have no immutable attempt binding')
        with binding_path.open('x') as stream:stream.write(json.dumps(binding,indent=2)+'\n')
    binding_sha=digest(binding_path)
    row_receipts={}
    for row in previous:
        shape=row['resolution'];path=root/f'{label}_{shape[0]}_{shape[1]}_{shape[2]}_row_binding.json'
        content=path.read_bytes()
        if json.loads(content)!=dict(attempt_binding_sha256=binding_sha,record=row):
            raise ValueError('retained row differs from its immutable evidence receipt')
        row_receipts[str(path)]=hashlib.sha256(content).hexdigest()
    def verify_attempt():
        verify_hashes(input_paths,input_hashes)
        if digest(binding_path)!=binding_sha:raise ValueError('immutable attempt binding changed')
        if any(digest(path)!=sha for path,sha in measured.items()):
            raise ValueError('measured CUDA build/dependencies changed during the attempt')
        if backend.dependency_images | loaded_kokkos_images()!=actual_images:
            raise ValueError('loaded puncture/runtime image changed during the attempt')
        if any(digest(path)!=sha for path,sha in row_receipts.items()):
            raise ValueError('retained row evidence receipt changed')
        for row in previous:
            artifacts=row.get('raw_artifact_sha256',{})
            expected={str((raw/f"{label}_{row['resolution'][0]}_{row['resolution'][2]}{suffix}.npz").resolve())
                      for suffix in ('','_collocation')}
            if set(artifacts)!=expected or any(digest(path)!=sha for path,sha in artifacts.items()):
                raise ValueError('retained row raw artifacts are missing or changed')
    verify_attempt()
    def factory(_backend,n,nphi):
        cfg=_backend.config();cfg.n[:]=grid
        spin=.99 if a.case.startswith('aligned') else 0.;speed=0. if spin else float(np.sqrt(.99))
        for h,sign in enumerate((1,-1)):
            cfg.hole[h]=Hole(.5,(sign*separation/2,0,0),(0,0,.25*spin),(-sign*speed,0,0))
            radius=.25*np.sqrt(1-spin*spin)*np.sqrt(1-speed*speed)
            cfg.inner_min[h]=.1*radius;cfg.inner_max[h]=.3*radius
        cfg.conformal_choice=0;cfg.inner_flatten=0;cfg.omega[:]=[1,1]
        cfg.attenuation_power=4;cfg.far_radius=0
        for field in ('max_newton','max_krylov','memory_limit_mib'):setattr(cfg,field,controls[field])
        cfg.tolerance=controls['outer_tolerance'];cfg.krylov_restart=controls['restart']
        return cfg
    result=solve_case(backend,factory,[(grid[0],grid[2])],label,horizon_scaled=True,adaptive_steps=True,
        previous_records=previous,execution='kokkos',solve_options=dict(krylov=controls['krylov'],linear_rtol=controls['linear_rtol']),
        output_report=report_path,raw_directory=raw)
    verify_attempt()
    record=result['records'][-1];cfg=factory(backend,grid[0],grid[2])
    row_path=root/f'{label}_{grid[0]}_{grid[1]}_{grid[2]}_row_binding.json'
    with row_path.open('x') as stream:
        stream.write(json.dumps(dict(attempt_binding_sha256=binding_sha,record=record),indent=2)+'\n')
    row_receipts[str(row_path)]=digest(row_path)
    artifacts=record['raw_artifact_sha256']
    if any(digest(path)!=sha for path,sha in artifacts.items()):raise ValueError('new row raw evidence changed')
    with np.load(raw/f'{label}_{grid[0]}_{grid[2]}.npz') as saved:values=saved['unknowns'].copy()
    if any(digest(path)!=sha for path,sha in artifacts.items()):raise ValueError('new row raw evidence changed while reading')
    exported=write_checkpoint(root/f'{label}_{grid[0]}_{grid[1]}_{grid[2]}.checkpoint',cfg,values,
                              backend.library_sha256(),'diagnostic',backend.parameterization())
    verify_attempt()
    result.update(stage='diagnostic_extreme_investigation',plan_sha256=input_hashes['plan'],
        performance_sha256=input_hashes['performance'],compiled_report_receipt=receipt,
        seed_controls_sha256=input_hashes['seed'],prerequisite_sha256=input_hashes,
        attempt_binding=dict(path=str(binding_path),sha256=binding_sha),prerequisites_passed=prerequisites_passed,
        library_dependency_images=backend.dependency_images,runtime_images=loaded_kokkos_images(),device=device,coordinate_separation=separation,
        portable_checkpoint=exported,acceptance_transferred=False,reference_reproduction=False,
        charge_quadrature_verified=False,coordinate_covariance_verified=False,horizon_enclosure_verified=False,
        binary_validation_complete=False)
    report=json.loads(report_path.read_text());report[label]=result
    verify_attempt()
    if any(digest(path)!=sha for path,sha in artifacts.items()):raise ValueError('new row raw evidence changed before saving')
    report_path.write_text(json.dumps(report,indent=2)+'\n')
    print(label,record['diagnostics'],exported,flush=True)
    return 0 if record['diagnostics']['status']==0 else 1


if __name__=='__main__':raise SystemExit(main())
