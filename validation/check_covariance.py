"""Solved binary rotation/translation and ADM-origin checks, one thread."""
import argparse,json,time,hashlib,os
from pathlib import Path
import numpy as np
from hispid import Backend,Config
from configs import moderate,as_dict
from physical import extrapolate
from run_validation import ROOT,REPORT,RAW,points,SOLVER_CONTROLS
from prolong import for_backend
from checkpoints import restore
from charge_checks import validate_refinements,qualify_refinement,measure_origin_quadrature,angular_changes

TENSORS=('gamma','Kij','conformal_metric','Atilde')

def transformed(config,Q,c):
    new=Config.from_buffer_copy(config)
    for h in new.hole:
        h.center[:]=Q@np.array(h.center)+c
        h.spin[:]=Q@np.array(h.spin);h.velocity[:]=Q@np.array(h.velocity)
    return new

def rotate_values(values,Q):
    out={k:v.copy() for k,v in values.items()}
    for k in TENSORS:out[k]=np.einsum('ik,nkl,jl->nij',Q,values[k].reshape(-1,3,3),Q).reshape(-1,9)
    out['correction'][:,1:]=values['correction'][:,1:]@Q.T
    return out

def errors(a,b):
    out={}
    for k in (*TENSORS,'psi','mean_curvature','correction'):
        diff=np.abs(a[k]-b[k]);scale=np.maximum(np.max(np.abs(b[k]),axis=1),1e-8) if b[k].ndim==2 else np.maximum(np.abs(b[k]),1e-8)
        out[k]=float(np.max(diff/scale[:,None])) if diff.ndim==2 else float(np.max(diff/scale))
    return out

def validate_source_sequence(records,levels):
    selected=[]
    for n,nphi in levels:
        candidates=[r for r in records if r['resolution'][0]==n and r['resolution'][2]==nphi]
        if len(candidates)!=1:raise ValueError('select a unique full grid for covariance')
        selected.append(candidates[0])
    shapes=[tuple(r['resolution']) for r in selected]
    if any(len(s)!=3 or list(s)!=r['config']['n'] for s,r in zip(shapes,selected)):
        raise ValueError('source record/configuration grid mismatch')
    if (len(set(shapes))!=len(shapes) or any(any(a>b for a,b in zip(old,new))
            for old,new in zip(shapes,shapes[1:]))):
        raise ValueError('full three-dimensional covariance grids must refine monotonically')
    data=[{k:v for k,v in r['config'].items() if k not in SOLVER_CONTROLS} for r in selected]
    if any(d!=data[0] for d in data[1:]) or len({r.get('horizon_scaled',False) for r in selected})!=1:
        raise ValueError('covariance truncation requires identical physical free data and sampling mode')
    return selected

def run(backend,levels,label='moderate',report_path=None,raw_directory=None,output_path=None,
        execution='reference',solve_options=None,charge_evidence=None,allow_diagnostic=False):
    source_path=REPORT if report_path is None else Path(report_path)
    source_bytes=source_path.read_bytes();source_sha=hashlib.sha256(source_bytes).hexdigest()
    saved=json.loads(source_bytes).get(label,{})
    options=dict(solve_options or {})
    from native_loader import validate_krylov,loaded_kokkos_images
    validate_krylov(options.get('krylov'),options.get('linear_rtol'))
    raw_root=RAW if raw_directory is None else Path(raw_directory)
    target=REPORT if output_path is None else Path(output_path)
    if report_path is not None and (output_path is None or target.resolve()==source_path.resolve()):
        raise ValueError('separate study results require a fresh separate covariance output')
    if output_path is not None and target.exists():raise FileExistsError('preserve prior covariance evidence')
    covariance_raw=raw_root if output_path is None else target.parent/'raw'/target.stem
    if (len(levels)<3 or len(set(levels))!=len(levels)
        or any(a>c or b>d for (a,b),(c,d) in zip(levels,levels[1:]))):
        raise ValueError('three distinct monotone covariance resolutions required')
    if (any((covariance_raw/f'covariance_{label}_{n}_{p}{suffix}.npz').exists()
            for n,p in levels for suffix in ('','_solve'))
        or (covariance_raw/f'translation_{label}_{levels[-1][0]}_{levels[-1][1]}_solve.npz').exists()):
        raise FileExistsError('covariance raw artifacts already exist; use a fresh output')
    selected_records=validate_source_sequence(saved.get('records',[]),levels)
    if not saved.get('passed') and not allow_diagnostic:raise ValueError('source physical convergence gate has not passed')
    explicit_study=bool(report_path is not None or raw_directory is not None or output_path is not None
        or execution!='reference' or options or charge_evidence or allow_diagnostic)
    source_stopping_verified=bool(not explicit_study or all(r.get('stopping_verified') is True for r in selected_records))
    charges_verified=bool(saved.get('charge_angular_verified') and not explicit_study)
    charge_bindings={};charge_per_grid={}
    if charge_evidence:
        charges_verified=True;shapes=[]
        for path in charge_evidence:
            payload=Path(path).read_bytes();charge=json.loads(payload);charge_bindings[str(Path(path).resolve())]=hashlib.sha256(payload).hexdigest()
            shapes.append(tuple(charge['resolution']))
            checked=qualify_refinement(charge['records'],charge['radii'],charge['qualification_bound'])
            charges_verified &= bool(charge.get('qualified_charge_checks') is True and charge['case']==label
                and charge['results_sha256']==source_sha and charge['library_sha256']==backend.library_sha256()
                and charge['unknown_parameterization_id']==backend.parameterization()
                and charge['collocation_maps']==backend.parameterization_maps() and checked['qualified_charge_checks'])
            charge_per_grid[tuple(charge['resolution'])]=charge
        selected=[tuple(r['resolution']) for r in selected_records]
        charges_verified &= sorted(shapes)==sorted(selected) and len(selected)==len(levels) and len(set(shapes))==len(shapes)
    if not charges_verified and not allow_diagnostic:raise ValueError('refined source charge checks have not passed')
    if not source_stopping_verified and not allow_diagnostic:raise ValueError('source stopping checks have not passed')
    prerequisites_passed=bool(saved.get('passed') and charges_verified and source_stopping_verified)
    sha=backend.library_sha256()
    if not saved.get('records') or any(r.get('library_sha256')!=sha for r in saved['records']):raise ValueError('moderate gate library SHA mismatch')
    frozen_images={str(backend.path):sha}|backend.dependency_images|loaded_kokkos_images()
    frozen_raw={}
    for record in saved['records']:
        if (record['resolution'][0],record['resolution'][2]) not in levels:continue
        artifacts=record.get('raw_artifact_sha256',{})
        if explicit_study and not artifacts:
            raise ValueError('explicit study covariance requires retained raw artifact hashes')
        if artifacts:
            for path,value in artifacts.items():frozen_raw[str(raw_root/Path(path).name)]=value
        else:
            path=raw_root/f"{record['case']}_{record['resolution'][0]}_{record['resolution'][2]}.npz"
            frozen_raw[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    saved_sha=hashlib.sha256(json.dumps(saved,sort_keys=True).encode()).hexdigest()
    def verify_inputs():
        current=source_path.read_bytes()
        if (output_path is not None and hashlib.sha256(current).hexdigest()!=source_sha
            or hashlib.sha256(json.dumps(json.loads(current).get(label,{}),sort_keys=True).encode()).hexdigest()!=saved_sha
            or backend.library_sha256()!=sha
            or any(hashlib.sha256(Path(path).read_bytes()).hexdigest()!=value for path,value in (frozen_images|charge_bindings|frozen_raw).items())):
            raise ValueError('bound covariance input or image changed')
    def save_progress(value):
        verify_inputs();target.parent.mkdir(parents=True,exist_ok=True)
        if output_path is None:
            report=json.loads(target.read_text());report['covariance']=value;target.write_text(json.dumps(report,indent=2)+'\n')
        else:target.write_text(json.dumps(value,indent=2)+'\n')
    def config_for(n,nphi):
        candidates=[r for r in saved['records'] if r['resolution'][0]==n and r['resolution'][2]==nphi]
        if len(candidates)!=1:raise ValueError('select a unique full grid for covariance')
        rec=candidates[0]
        return restore(backend,rec,raw_root)
    def stop_verified(diagnostic,cfg):
        return bool(diagnostic['status']==0 and diagnostic['converged']
                    and np.isfinite(diagnostic['scaled_linf']).all() and max(diagnostic['scaled_linf'])<=cfg.tolerance)
    axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis);angle=.73
    W=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    Q=np.eye(3)+np.sin(angle)*W+(1-np.cos(angle))*(W@W);offset=np.array([.7,-.2,.4])
    records=[];attempts=[];base_values=[];previous_values=None;previous_shape=None
    def progress():
        return dict(source_case=label,records=records,attempts=attempts,passed=False,
            source_results_sha256=source_sha,prerequisites_passed=prerequisites_passed,
            source_stopping_verified=source_stopping_verified,explicit_study=explicit_study,
            bound_images=frozen_images,charge_evidence_sha256=charge_bindings,
            source_raw_artifact_sha256=frozen_raw,unknown_parameterization_id=backend.parameterization(),
            collocation_maps=backend.parameterization_maps(),acceptance_inherited=False,binary_validation_complete=False)
    def retain_solve(solution,attempt,config,unknowns,path):
        verify_inputs();path.parent.mkdir(parents=True,exist_ok=True)
        if path.exists():raise FileExistsError('preserve completed solve coefficients')
        np.savez_compressed(path,unknowns=unknowns)
        attempt.update(config=as_dict(config),library_sha256=sha,
            solve_raw_artifact=dict(path=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        if explicit_study:
            from execution import name,concurrency,device_description
            attempt.update(compiled_execution=name(backend.lib),execution_concurrency=concurrency(backend.lib),
                device=device_description(backend.lib),resolved_solve_options=solution.resolved_options,
                linear_history=solution.linear_history(),work_statistics=solution.work_statistics())
        save_progress(progress())
    save_progress(progress())
    for n,nphi in levels:
        verify_inputs();config,base_unknowns=config_for(n,nphi)
        base_record=next(r for r in saved['records'] if r['resolution']==list(config.n))
        attempt=dict(kind='rotation',resolution=list(config.n),stage='sampling_base',completed=False)
        attempts.append(attempt);save_progress(progress())
        x,near,bulk=points(config,base_record.get('horizon_scaled',False));x=x[:near+bulk]
        start=time.monotonic()
        with backend.create_sampler(config) as s:
            s.set_unknowns(base_unknowns)
            original=s.sample(x);base_values.append(original)
        prime=transformed(config,Q,offset)
        attempt['stage']='solving';save_progress(progress())
        with backend.create(prime,execution=execution) as t:
            if previous_values is not None:t.set_unknowns(for_backend(backend,previous_values,previous_shape,list(prime.n)))
            diagnostic=t.solve(**options)
            attempt.update(stage='measuring_fields_and_charges',diagnostics=diagnostic,
                stopping_verified=stop_verified(diagnostic,prime))
            previous_values=t.unknowns();previous_shape=list(prime.n)
            retain_solve(t,attempt,prime,previous_values,covariance_raw/f'covariance_{label}_{n}_{nphi}_solve.npz')
            rotated=t.sample(x@Q.T+offset)
            err=errors(rotated,rotate_values(original,Q))
            source_charge=charge_per_grid.get(tuple(config.n))
            refined_charge=dict(qualified_charge_checks=False,records=[])
            if source_charge is not None:
                radii=source_charge['radii'];quadratures=[(r['ntheta'],r['nphi']) for r in source_charge['records']]
                bound=source_charge['qualification_bound'];validate_refinements(radii,quadratures,bound)
                refined_charge.update(radii=radii,quadratures=quadratures,qualification_bound=bound)
                for nt,np_ in quadratures:
                    item=measure_origin_quadrature(t,radii,nt,np_,offset)
                    angular_changes(item,refined_charge['records'][-1] if refined_charge['records'] else None,
                        prefixes=('', 'native_', 'global_', 'global_native_'),global_origin=True)
                    refined_charge['records'].append(item);attempt['charge_refinement']=refined_charge
                    save_progress(progress())
                refined_charge.update(qualify_refinement(refined_charge['records'],radii,bound,
                    global_origin=True,sampled_rotation=False))
                qrot=np.asarray(refined_charge['records'][-1]['extrapolated_EPJ'])
                qglobal=np.asarray(refined_charge['records'][-1]['global_extrapolated_EPJ'])
                base_charge=np.asarray(source_charge['records'][-1]['independent_extrapolated_EPJ'])
            else:
                radii=[100.,200.,400.]
                qrot=extrapolate(radii,[t.charges(R,center=offset,ntheta=12,nphi=24) for R in radii])
                qglobal=extrapolate(radii,[t.charges(R,ntheta=12,nphi=24) for R in radii])
                base_charge=np.asarray(base_record['charges_extrapolated'])
            execution_evidence={key:attempt[key] for key in ('compiled_execution','execution_concurrency','device',
                'resolved_solve_options','linear_history','work_statistics')} if explicit_study else {}
        expected=np.r_[base_charge[0],Q@base_charge[1:4],Q@base_charge[4:]]
        origin_expected=expected.copy();origin_expected[4:]+=np.cross(offset,expected[1:4])
        verify_inputs();covariance_raw.mkdir(parents=True,exist_ok=True)
        artifact=covariance_raw/f'covariance_{label}_{n}_{nphi}.npz'
        np.savez_compressed(artifact,points=x,rotated_unknowns=previous_values,
            **{'base_'+k:v for k,v in original.items()},**{'rotated_'+k:v for k,v in rotated.items()})
        rec=dict(resolution=list(prime.n),library_sha256=backend.library_sha256(),rotation=Q.tolist(),offset=offset.tolist(),config=as_dict(prime),diagnostics=diagnostic,
                 relative_errors=err,rotated_EPJ=qrot.tolist(),global_origin_EPJ=qglobal.tolist(),expected_EPJ=expected.tolist(),
                 charge_error=np.abs(qrot-expected).tolist(),origin_charge_error=np.abs(qglobal-origin_expected).tolist(),seconds=time.monotonic()-start,
                 execution=execution,requested_solve_options=options,stopping_verified=stop_verified(diagnostic,prime),
                 charge_refinement=refined_charge,
                 raw_artifact=dict(path=str(artifact.resolve()),sha256=hashlib.sha256(artifact.read_bytes()).hexdigest()))
        rec.update(execution_evidence)
        records.append(rec);attempt.update(stage='completed',completed=True);print(json.dumps(rec),flush=True)
        save_progress(progress())
    truncation=errors(base_values[-2],base_values[-1]);last=records[-1]
    field_pass=all(last['relative_errors'][k]<=max(1e-9,5*truncation[k]) for k in truncation)
    charge_accuracy_verified=bool(all(r['charge_refinement']['qualified_charge_checks'] for r in records))
    charge_pass=max(last['charge_error']+last['origin_charge_error'])<.005 and (charge_accuracy_verified or not explicit_study)
    # Translation alone leaves the canonical-frame coefficients invariant.
    verify_inputs();n,nphi=levels[-1];cfg,base_unknowns=config_for(n,nphi)
    best=next(r for r in saved['records'] if r['resolution']==list(cfg.n))
    attempt=dict(kind='translation',resolution=list(cfg.n),stage='solving',completed=False)
    attempts.append(attempt);save_progress(progress())
    x,near,bulk=points(cfg,best.get('horizon_scaled',False));x=x[:near+bulk]
    with backend.create(transformed(cfg,np.eye(3),offset),execution=execution) as t:
        t.set_unknowns(base_unknowns);d=t.solve(**options)
        translation_cfg=transformed(cfg,np.eye(3),offset)
        attempt.update(stage='measuring_fields',diagnostics=d,stopping_verified=stop_verified(d,cfg))
        retain_solve(t,attempt,translation_cfg,t.unknowns(),covariance_raw/f'translation_{label}_{n}_{nphi}_solve.npz')
        translation_errors=errors(t.sample(x+offset),base_values[-1])
    attempt.update(stage='completed',completed=True);save_progress(progress())
    translation_pass=max(translation_errors.values())<1e-9 and (stop_verified(d,cfg) if explicit_study else d['status']==0)
    verify_inputs()
    for rec in records:
        if hashlib.sha256(Path(rec['raw_artifact']['path']).read_bytes()).hexdigest()!=rec['raw_artifact']['sha256']:
            raise ValueError('covariance raw evidence changed')
    for attempt in attempts:
        artifact=attempt['solve_raw_artifact']
        if hashlib.sha256(Path(artifact['path']).read_bytes()).hexdigest()!=artifact['sha256']:
            raise ValueError('retained fresh solve coefficients changed')
    measured_pass=bool(field_pass and charge_pass and translation_pass and all(
        r['stopping_verified'] if explicit_study else r['diagnostics']['status']==0 for r in records))
    result=dict(source_case=label,records=records,attempts=attempts,finest_base_truncation=truncation,translation_errors=translation_errors,translation_diagnostics=d,
                field_pass=bool(field_pass),charge_pass=bool(charge_pass),translation_pass=bool(translation_pass),
                passed=bool(measured_pass and prerequisites_passed),measured_covariance_checks_passed=measured_pass,
                prerequisites_passed=prerequisites_passed,source_results_sha256=source_sha,bound_images=frozen_images,
                source_stopping_verified=source_stopping_verified,explicit_study=explicit_study,
                charge_evidence_sha256=charge_bindings,source_raw_artifact_sha256=frozen_raw,
                charge_accuracy_verified=charge_accuracy_verified,
                acceptance_inherited=False,binary_validation_complete=False)
    if output_path is not None:save_progress(result)
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--library',required=True);p.add_argument('--levels',default='24:12,40:20,56:28');p.add_argument('--case',default='moderate')
    p.add_argument('--results');p.add_argument('--raw-directory');p.add_argument('--output')
    p.add_argument('--charge-evidence',action='append',help='repeat once per source grid with qualified independent/angular/radial checks')
    p.add_argument('--allow-diagnostic',action='store_true');p.add_argument('--execution',choices=('reference','kokkos'),default='reference')
    p.add_argument('--threads',type=int,default=1);p.add_argument('--krylov',choices=('gmres','bicgstab'));p.add_argument('--linear-rtol',type=float);a=p.parse_args()
    if not 1<=a.threads<=16:raise ValueError('host concurrency must be1..16')
    os.environ.update(OMP_NUM_THREADS=str(a.threads),OPENBLAS_NUM_THREADS='1')
    b=Backend(a.library)
    from execution import select,name,device_description
    select(b.lib,a.execution,a.threads)
    if a.execution=='kokkos' and name(b.lib)=='Cuda':
        device=device_description(b.lib)
        if not device or device['visible_count']!=1 or device['visible_ordinal']!=0:raise ValueError('one allocated GPU required')
    options={key:value for key,value in dict(krylov=a.krylov,linear_rtol=a.linear_rtol).items() if value is not None}
    r=run(b,[tuple(map(int,s.split(':'))) for s in a.levels.split(',')],a.case,
        report_path=a.results,raw_directory=a.raw_directory,output_path=a.output,execution=a.execution,
        solve_options=options,charge_evidence=a.charge_evidence,allow_diagnostic=a.allow_diagnostic)
    if a.output is None:
        report=json.loads(REPORT.read_text());report['covariance']=r;REPORT.write_text(json.dumps(report,indent=2)+'\n')
    return 0 if r['passed'] else 1

if __name__=='__main__':raise SystemExit(main())
