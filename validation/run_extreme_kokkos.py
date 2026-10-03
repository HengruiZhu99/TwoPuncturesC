"""One retained grid of a separate aligned-spin or head-on boost investigation.

Normally requires the completed performance matrix and a compiled-report
receipt. An explicitly recorded human override can waive that sequencing.
Uses the common physical validation machinery with
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


def human_override_authorized(override,performance_sha):
    if override is None:return False
    if (override.get('schema')!='hispid_human_override_v1'
        or override.get('authority')!='explicit_user_request'
        or override.get('scope')!='stop_performance_and_proceed_to_extreme_investigations'
        or override.get('performance_sha256')!=performance_sha
        or not override.get('instruction') or override.get('scientific_acceptance_waived') is not False):
        raise ValueError('human override must record the explicit scope and preserved performance snapshot')
    return True


def prerequisites(performance,receipt,hashes,human_override=None):
    grids=[[40,80,16],[80,160,16],[128,256,28]]
    variants={'reference','serial','openmp1','openmp2','openmp4','openmp8','openmp16','cuda'}
    systems={'hispid','by'};methods={'gmres','bicgstab'}
    binding=performance.get('binding',{})
    ids=[v['id'] for v in performance.get('manifest',{}).get('variants',[])]
    expected={f'{"_".join(map(str,grid))}_{variant}_{system}_{method}_{repeat}'
        for grid in grids for variant in variants for system in systems for method in methods for repeat in range(3)}
    records=performance.get('records',{});failures=performance.get('failures',{})
    overridden=human_override_authorized(human_override,hashes['performance'])
    if ((not overridden and (performance.get('declared_performance_completed') is not True
        or performance.get('expected_workers')!=288 or performance.get('completed_workers')!=288
        or set(records)|set(failures)!=expected))
        or binding.get('grids')!=grids or binding.get('repeats')!=3
        or len(binding.get('systems',[]))!=2 or set(binding.get('systems',[]))!=systems
        or len(binding.get('methods',[]))!=2 or set(binding.get('methods',[]))!=methods
        or len(ids)!=8 or set(ids)!=variants or set(records)&set(failures)
        or (set(records)|set(failures))-expected):
        raise ValueError('complete declared matrix required unless explicitly overridden; retained identities must remain valid')
    for label,row in records.items():
        actual=f'{"_".join(map(str,row.get("grid",[])))}_{row.get("variant")}_{row.get("mode")}_{row.get("krylov")}_{row.get("repeat")}'
        if actual!=label:raise ValueError('performance record identity differs from its matrix key')
    def valid_sha(value):
        return isinstance(value,str) and len(value)==64 and all(c in '0123456789abcdef' for c in value)
    for rows,keys in ((records,('state','log','worker')),(failures,('log',))):
        for row in rows.values():
            if any(not isinstance(row.get(key),str) or not row[key]
                   or not valid_sha(row.get(key+'_sha256')) for key in keys):
                raise ValueError('complete retained performance artifact witnesses required')
    inputs=performance.get('input_sha256',{})
    required_inputs={f'input_{"_".join(map(str,grid))}.json' for grid in grids}
    if (len(inputs)!=3 or any(not isinstance(path,str) or not path for path in inputs)
        or {Path(path).name for path in inputs}!=required_inputs
        or any(not valid_sha(sha) for sha in inputs.values())):
        raise ValueError('all three retained performance input witnesses required')
    if not overridden and (receipt.get('compilation_confirmed') is not True
        or receipt.get('compiler')!='mcp__codex_app__compile_latex_document'
        or receipt.get('performance_sha256')!=hashes['performance']
        or receipt.get('report_sha256')!=hashes['report']):
        raise ValueError('matching successful native LaTeX compilation receipt required')
    return performance,receipt


def completed_seed_controls(seed,geometry='host'):
    """Numerical failures may be diagnostic; unfinished controls are not ready."""
    cases=seed.get('cases',[])
    if (seed.get('completed') is not True or len(cases)!=2
        or {row.get('case') for row in cases}!={'spin99','gamma10'}
        or any(row.get('completed') is not True for row in cases)):
        raise ValueError('completed fresh chi=.99/Gamma=10 seed controls required')
    if geometry not in ('host','execution'):raise ValueError('invalid seed geometry selection')
    if geometry=='execution':
        device=seed.get('device') or {}
        if (seed.get('seed_execution')!='kokkos' or seed.get('compiled_execution')!='Cuda'
            or seed.get('seed_scalar_digits')!=53 or seed.get('exact_horizon_geometry')!='host'
            or device.get('visible_count')!=1 or device.get('visible_ordinal')!=0 or not device.get('uuid')):
            raise ValueError('execution geometry requires fresh device-field seed controls and separately labelled host horizon gradients')
    return seed


def verify_seed_images(seed,measured,producer,geometry='host'):
    """Bind seed controls to the measured producer, puncture and runtime DSOs."""
    bound=seed.get('bound_artifacts_sha256',{})
    def native_path(path):
        return Path(path).name.startswith(('libHiSpID','libTwoPunctures','libkokkos'))
    images=seed.get('native_images') or {p:s for p,s in bound.items() if native_path(p)}
    # Historical host controls without an image witness remain compatible.
    # Explicit execution geometry always needs the complete native witnesses.
    if not images and geometry=='host':return
    required={str(producer)}|{p for p in measured if Path(p).name.startswith(('libTwoPunctures','libkokkos'))}
    dependencies={**seed.get('library_dependency_images',{}),**seed.get('runtime_images',{})}
    if (not required.issubset(images) or not any('libkokkos' in Path(p).name for p in images)
        or any(measured.get(p)!=sha or bound.get(p)!=sha for p,sha in images.items())
        or any(images.get(p)!=sha for p,sha in dependencies.items())):
        raise ValueError('seed producer, puncture or Kokkos dependency images differ from measured evidence')


def source_floor_passed(floor,library_sha,grids,separation,mass=.5,memory_mib=32768,geometry='host',target_case=None):
    """Bind exact isolated controls to the binary's charts and declared norm.

    Structural mismatches are rejected; retained numerical failures return
    false and may only accompany an explicitly diagnostic investigation.
    """
    if (geometry not in ('host','execution') or floor.get('geometry','host')!=geometry
        or floor.get('schema')!='hispid_source_floor_v3' or floor.get('extreme_controls') is not True
        or floor.get('library_sha256')!=library_sha or floor.get('residual_scaling')!='sin3_alpha_beta'
        or floor.get('execution')!='kokkos' or floor.get('compiled_execution')!='Cuda'
        or floor.get('seed_mass')!=mass or floor.get('coordinate_separation')!=separation
        or floor.get('weighted_far_source_bound')!=1e-14 or floor.get('far_radius_minimum')!=100.
        or floor.get('exact_correction')!=0.):
        raise ValueError('source-floor controls differ from the target charts, norm or producer')
    from check_far_source_floor import isolated_controls
    controls={label:(spin*mass**2,velocity) for label,spin,velocity in isolated_controls(True)}
    if target_case is not None:
        if target_case not in ('spin99','gamma10') or floor.get('target_case') not in (None,target_case):
            raise ValueError('source-floor target case differs from the binary')
        controls={target_case:controls[target_case]}
    rows=floor.get('records',[])
    keys=[(r['case'],r['active_hole'],tuple(r['resolution'])) for r in rows]
    required={(label,active,tuple(grid)) for label in controls for active in (0,1) for grid in grids}
    if len(keys)!=len(set(keys)) or not required.issubset(keys):
        raise ValueError('complete distinct exact-seed controls required at every target grid')
    passed=floor.get('passed') is True
    for row in rows:
        if (row['case'],row['active_hole'],tuple(row['resolution'])) not in required:continue
        setup=row.get('setup_statistics')
        if setup and setup.get('geometry_execution')!=int(geometry=='execution'):
            raise ValueError('source-floor setup witness contradicts the selected geometry')
        if geometry=='execution' and (not setup or setup.get('geometry_execution')!=1 or setup.get('scalar_digits')!=53):
            raise ValueError('source-floor execution geometry witness differs from the binary')
        cfg=row['config'];spin,velocity=controls[row['case']]
        active=row['active_hole'];sign=1 if active==0 else -1
        hole=cfg['hole'][active];inactive=cfg['hole'][1-active]
        if (cfg['n']!=row['resolution'] or cfg['memory_limit_mib']!=memory_mib
            or hole['mass']!=mass or hole['center']!=[sign*separation/2,0,0]
            or not np.array_equal(hole['spin'],spin) or not np.array_equal(hole['velocity'],-sign*velocity)
            or inactive['mass']!=0 or inactive['center']!=[-sign*separation/2,0,0]
            or cfg['conformal_choice']!=0 or cfg['inner_flatten']!=0
            or cfg['inner_min']!=[0,0] or cfg['inner_max']!=[0,0]
            or cfg['omega']!=[0,0] or cfg['far_radius']!=0):
            raise ValueError('source-floor geometry or attenuation differs from exact isolated controls')
        values=np.asarray(row.get('weighted_far_source_linf',[]))
        physical_values=np.asarray(row.get('physical_equivalent_far_linf',[]))
        passed &= bool(row.get('completed') is True and row.get('passed') is True and values.shape==(4,)
            and np.isfinite(values).all() and np.min(values)>=0 and np.max(values)<1e-14
            and physical_values.shape==(4,) and np.isfinite(physical_values).all() and np.min(physical_values)>=0
            and row.get('far_node_count',0)>0 and row.get('far_radius_range',[0])[0]>=100.)
    return bool(passed)


def verify_source_floor_arrays(floor):
    """Recompute admission summaries from their immutable full numerical arrays."""
    bindings={}
    for row in floor['records']:
        artifact=row['raw_artifact'];path=artifact['path'];sha=artifact['sha256']
        if path in bindings:raise ValueError('each source-floor control needs a distinct retained artifact')
        if digest(path)!=sha:raise ValueError('source-floor array hash mismatch')
        with np.load(path) as data:
            xyz=data['xyz'];weighted=data['weighted'];unknowns=data['unknowns'];saved_mask=data['far_mask']
            physical=data['physical_equivalent']
            count=int(np.prod(row['resolution']))
            if (xyz.shape!=(count,3) or weighted.shape!=(count,4) or physical.shape!=(count,4)
                or unknowns.size!=4*count or not np.isfinite(xyz).all() or not np.isfinite(unknowns).all()
                or not np.all(unknowns==0) or saved_mask.shape!=(count,)):
                raise ValueError('source-floor array dimensions or zero-correction witness differ')
            radius=np.linalg.norm(xyz,axis=1);mask=radius>=floor['far_radius_minimum']
            if not mask.any() or not np.array_equal(saved_mask,mask):raise ValueError('source-floor far mask differs')
            values=np.max(abs(weighted[mask]),axis=0);pvalues=np.max(abs(physical[mask]),axis=0)
            limits=[float(radius[mask].min()),float(radius[mask].max())]
            if (int(mask.sum())!=row['far_node_count'] or not np.array_equal(values,row['weighted_far_source_linf'],equal_nan=True)
                or not np.array_equal(pvalues,row['physical_equivalent_far_linf'],equal_nan=True)
                or not np.allclose(limits,row['far_radius_range'],rtol=4*np.finfo(float).eps,atol=0)):
                raise ValueError('source-floor retained arrays differ from reported summaries')
        if digest(path)!=sha:raise ValueError('source-floor arrays changed while decoding')
        bindings[path]=sha
    return bindings


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library',required=True);p.add_argument('--performance-results',required=True)
    p.add_argument('--performance-artifact-root',required=True,help='original benchmark source root containing retained state/log/worker/input files')
    p.add_argument('--compiled-report-receipt');p.add_argument('--compiled-report')
    p.add_argument('--human-override',help='Recorded explicit user instruction to stop performance work and proceed; never waives scientific acceptance')
    p.add_argument('--seed-controls',required=True);p.add_argument('--source-floor-controls',required=True)
    p.add_argument('--plan',default=str(Path(__file__).with_name('extreme_kokkos_plan.json')))
    p.add_argument('--case',choices=('aligned_spin99_kokkos','headon_gamma10_kokkos'),required=True)
    p.add_argument('--grid-index',type=int,required=True)
    p.add_argument('--separation',type=float,help='retain a separately labeled separation-calibration case')
    p.add_argument('--output-directory',required=True);p.add_argument('--threads',type=int,default=16)
    p.add_argument('--geometry',choices=('host','execution'),default='host')
    p.add_argument('--allow-diagnostic-investigation',action='store_true',help='retain failed prerequisite gates while measuring unqualified data')
    a=p.parse_args()
    if not a.human_override and not (a.compiled_report_receipt and a.compiled_report):
        p.error('compiled report and receipt required without a human override')
    if bool(a.compiled_report_receipt)!=bool(a.compiled_report):p.error('provide both report and receipt, or neither under the override')
    input_paths=dict(performance=a.performance_results,plan=a.plan,seed=a.seed_controls,floor=a.source_floor_controls)
    if a.compiled_report:input_paths.update(receipt=a.compiled_report_receipt,report=a.compiled_report)
    if a.human_override:input_paths['human_override']=a.human_override
    inputs,input_hashes=frozen_inputs(input_paths)
    performance,receipt=prerequisites(inputs['performance'],inputs.get('receipt',{}),input_hashes,inputs.get('human_override'))
    from benchmark_kokkos import verify_artifacts
    performance_artifact_root=Path(a.performance_artifact_root).resolve(strict=True)
    verify_artifacts(performance,performance_artifact_root)
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
    seed=completed_seed_controls(inputs['seed'],geometry=a.geometry)
    verify_seed_images(seed,measured,library,geometry=a.geometry)
    if seed.get('library_sha256')!=digest(library) or {c['case'] for c in seed.get('cases',[])}!={'spin99','gamma10'}:
        raise ValueError('fresh bound chi=.99/Gamma=10 seed controls required')
    seed_cases={c['case']:c for c in seed['cases']}
    if seed_cases['spin99']['seed_rest_chi']!=.99 or abs(seed_cases['gamma10']['input_lorentz_factor']/10-1)>1e-13:
        raise ValueError('seed controls do not represent the requested targets')
    initial_separation=12. if a.case.startswith('aligned') else 25.
    separation=a.separation if a.separation is not None else initial_separation
    floor=inputs['floor']
    floor_passed=source_floor_passed(floor,digest(library),grids,separation,memory_mib=controls['memory_limit_mib'],geometry=a.geometry,
                                    target_case='spin99' if a.case.startswith('aligned') else 'gamma10')
    floor_artifacts=verify_source_floor_arrays(floor)
    if (not floor.get('bound_images') or any(measured.get(path)!=sha for path,sha in floor['bound_images'].items())
        or any(digest(path)!=sha for path,sha in floor_artifacts.items())):
        raise ValueError('source-floor images or retained arrays differ from measured evidence')
    prerequisites_passed=bool(seed.get('passed') and performance.get('all_stopping_checks_passed')
                              and performance.get('strict_state_comparisons_passed') and floor_passed)
    if not prerequisites_passed and not a.allow_diagnostic_investigation:
        raise ValueError('failed prerequisites require an explicitly diagnostic investigation')
    os.environ.update(OMP_NUM_THREADS=str(a.threads),OPENBLAS_NUM_THREADS='1',OMP_PROC_BIND='close',OMP_PLACES='cores')
    backend=Backend(str(library));select(backend.lib,'kokkos',a.threads)
    if (backend.parameterization_maps()!=dict(radial_stretch=plan['common_free_data']['maps'][0],angular_stretch=plan['common_free_data']['maps'][1])
        or backend.residual_scaling()!='sin3_alpha_beta' or plan['solve_controls']['row_power']!=3):
        raise ValueError('producer basis/maps/row scaling differs from the declared plan')
    if floor.get('unknown_parameterization_id')!=backend.parameterization() or floor.get('collocation_maps')!=backend.parameterization_maps():
        raise ValueError('source-floor continuous basis/maps differ from producer')
    actual_images=backend.dependency_images | loaded_kokkos_images()
    if any(measured.get(path)!=sha for path,sha in actual_images.items()):
        raise ValueError('loaded puncture/runtime image differs from the measured CUDA build')
    device=device_description(backend.lib)
    if name(backend.lib)!='Cuda' or concurrency(backend.lib)!=a.threads or not device or device['visible_count']!=1 or device['visible_ordinal']!=0:
        raise ValueError('actual one-GPU CUDA execution and requested host concurrency required')
    label=a.case+(('_d'+format(separation,'.17g').replace('.','p')) if a.separation is not None else '')
    fourier=a.grid_index==len(case['grids'])
    if fourier:label+='_fourier_control'
    root=Path(a.output_directory).resolve();root.mkdir(parents=True,exist_ok=True)
    report_path=root/'results.json';raw=root/'raw';raw.mkdir(exist_ok=True)
    report=json.loads(report_path.read_text()) if report_path.exists() else {}
    previous=report.get(label,{}).get('records',[])
    grid=grids[a.grid_index]
    if any(r['resolution']==grid for r in previous):raise ValueError('grid already retained; use a separate case for a new attempt')
    if any((raw/f'{label}_{grid[0]}_{grid[2]}{suffix}.npz').exists() for suffix in ('','_collocation','_solve')):
        raise FileExistsError('unrecorded attempt artifacts exist; preserve them and use a fresh output location')
    if not fourier and a.grid_index!=len(previous):raise ValueError('retain the declared radial sequence in order')
    if [r['resolution'] for r in previous]!=case['grids'][:len(previous)]:
        raise ValueError('retained resolutions differ from the declared grid prefix')
    binding=dict(schema='hispid_extreme_attempt_binding_v1',case=a.case,label=label,
        prerequisite_sha256=input_hashes,measured_images=measured,loaded_images=actual_images,
        producer=str(library),coordinate_separation=separation,host_threads=a.threads,geometry=a.geometry,
        human_override=inputs.get('human_override'))
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
        if any(digest(path)!=sha for path,sha in floor_artifacts.items()):
            raise ValueError('source-floor arrays changed during the attempt')
        if backend.dependency_images | loaded_kokkos_images()!=actual_images:
            raise ValueError('loaded puncture/runtime image changed during the attempt')
        if any(digest(path)!=sha for path,sha in row_receipts.items()):
            raise ValueError('retained row evidence receipt changed')
        for row in previous:
            artifacts=row.get('raw_artifact_sha256',{})
            expected={str((raw/f"{label}_{row['resolution'][0]}_{row['resolution'][2]}{suffix}.npz").resolve())
                      for suffix in ('','_collocation','_solve')}
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
        previous_records=previous,execution='kokkos',geometry=a.geometry,solve_options=dict(krylov=controls['krylov'],linear_rtol=controls['linear_rtol']),
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
        human_override=inputs.get('human_override'),performance_sequence_overridden=bool(a.human_override),
        performance_artifact_root=str(performance_artifact_root),
        seed_controls_sha256=input_hashes['seed'],prerequisite_sha256=input_hashes,
        source_floor_controls_sha256=input_hashes['floor'],source_floor_passed=floor_passed,
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
