"""Bound adaptive FD calibration of a separate study, without acceptance transfer.

The saved binary, exact isolated seeds and analytic Brill--Lindquist vacuum
data use identical points and pointwise steps. Raw convergence failures stay
unchanged; this replay cannot qualify a binary or waive an exterior stencil.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from hispid import Backend
from checkpoints import restore
from check_covariance import validate_source_sequence
from run_validation import points
from physical import constraints,norms
from native_loader import loaded_kokkos_images


def sampling_plan(config,record,region,factors):
    if (len(factors)<3 or not np.isfinite(factors).all() or min(factors)<=0
        or not np.all(np.diff(factors)<0)):
        raise ValueError('at least three finite positive decreasing step factors required')
    xyz,near,bulk=points(config,record.get('horizon_scaled',False))
    steps=np.asarray(record['verifier_steps'])
    if steps.shape!=(len(xyz),) or not np.isfinite(steps).all() or np.min(steps)<=0:
        raise ValueError('saved pointwise verifier steps do not match the generated points')
    roles=np.r_[np.zeros(near,dtype=int),np.ones(bulk,dtype=int),np.full(len(xyz)-near-bulk,2,dtype=int)]
    selected={'all':np.ones(len(xyz),dtype=bool),'exterior':roles<2,'bulk':roles==1}[region]
    return xyz[selected],roles[selected],[factor*steps[selected] for factor in factors]


def verify_saved_plan(config,record,raw_directory):
    """Bind the generated full plan to the original, hash-retained sampling."""
    n,_,nphi=record['resolution']
    path=Path(raw_directory)/f"{record['case']}_{n}_{nphi}.npz"
    bound=[sha for file,sha in record.get('raw_artifact_sha256',{}).items()
           if Path(file).name==path.name]
    if len(bound)!=1 or hashlib.sha256(path.read_bytes()).hexdigest()!=bound[0]:
        raise ValueError('uniquely hash-bound verifier sampling payload required')
    with np.load(path,allow_pickle=False) as saved:
        xyz=saved['points'].copy();steps=saved['verifier_steps'].copy()
    if hashlib.sha256(path.read_bytes()).hexdigest()!=bound[0]:
        raise ValueError('verifier sampling payload changed during decode')
    generated,_,_=points(config,record.get('horizon_scaled',False))
    if not (np.array_equal(xyz,generated) and np.array_equal(steps,record['verifier_steps'])):
        raise ValueError('saved source points or steps disagree with the calibration plan')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('library','results','raw-directory','case','output'):p.add_argument('--'+name,required=True)
    p.add_argument('--region',choices=('all','exterior','bulk'),default='all')
    p.add_argument('--step-factors',default='2,1,.5,.25,.125');a=p.parse_args()
    target=Path(a.output)
    if target.exists():raise FileExistsError('preserve previous verifier calibration')
    source=Path(a.results);source_bytes=source.read_bytes();source_sha=hashlib.sha256(source_bytes).hexdigest()
    case=json.loads(source_bytes)[a.case];records=case['records']
    if not records:raise ValueError('a retained source iterate is required')
    validate_source_sequence(records,[(r['resolution'][0],r['resolution'][2]) for r in records])
    backend=Backend(a.library);cfg,_=restore(backend,records[-1],a.raw_directory)
    verify_saved_plan(cfg,records[-1],a.raw_directory)
    factors=list(map(float,a.step_factors.split(',')))
    xyz,roles,steps=sampling_plan(cfg,records[-1],a.region,factors)
    frozen={str(backend.path):backend.library_sha256()}|backend.dependency_images|loaded_kokkos_images()
    for record in records:
        if (record['library_sha256']!=backend.library_sha256() or not record.get('raw_artifact_sha256')
            or record.get('unknown_parameterization_id')!=backend.parameterization()
            or record.get('unknown_parameterization')!=backend.parameterization_description()
            or record.get('collocation_maps')!=backend.parameterization_maps()):
            raise ValueError('every study grid requires matching library and retained raw hashes')
        for file,sha in record['raw_artifact_sha256'].items():frozen[str(Path(a.raw_directory)/Path(file).name)]=sha
    output=dict(case=a.case,results_sha256=source_sha,source_results=str(source.resolve()),
        library_sha256=backend.library_sha256(),bound_files_sha256=frozen,
        unknown_parameterization_id=backend.parameterization(),collocation_maps=backend.parameterization_maps(),
        region=a.region,points=xyz.tolist(),point_roles=roles.tolist(),step_factors=factors,
        pointwise_steps=[step.tolist() for step in steps],controls=[],binary_grids=[],
        completed=False,acceptance=False,changes_original_gate=False,acceptance_inherited=False,
        note='A noise calibration cannot relabel the source raw convergence flags or pass a binary.')
    def verify():
        if hashlib.sha256(source.read_bytes()).hexdigest()!=source_sha or any(
            hashlib.sha256(Path(file).read_bytes()).hexdigest()!=sha for file,sha in frozen.items()):
            raise ValueError('bound verifier source, array or image changed')
    def save():
        verify();target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(output,indent=2)+'\n')
    def sequence(sample,item):
        item['sequence']=[];save()
        for factor,step in zip(factors,steps):
            verify();result=constraints(sample,xyz,step)
            item['sequence'].append(dict(step_factor=factor,norms=norms(result),
                point_regions={label:norms(result,roles==code) for code,label in enumerate(('near','bulk','modified_sample_points'))},
                stencil_g_equals_one=norms(result,result['stencil_attenuation_all_one']),
                stencil_touches_modified=norms(result,~result['stencil_attenuation_all_one']),
                exterior_stencil_verified=bool(np.all(result['stencil_attenuation_all_one'][roles<2])),
                min_metric_eigenvalue=float(np.min(result['min_metric_eigenvalue'])),
                H=result['H'].tolist(),M=result['M'].tolist(),attenuation=result['attenuation'].tolist(),
                stencil_attenuation_min=result['stencil_attenuation_min'].tolist()))
            save()
    def brill_lindquist(x):
        x=np.asarray(x);psi=np.ones(len(x))
        for hole in cfg.hole:
            if hole.mass>0:psi+=hole.mass/(2*np.linalg.norm(x-np.array(hole.center),axis=1))
        return dict(gamma=(psi[:,None,None]**4*np.eye(3)).reshape(-1,9),Kij=np.zeros((len(x),9)),attenuation=np.ones(len(x)))
    save();item=dict(case='exact analytic Brill-Lindquist');output['controls'].append(item);sequence(brill_lindquist,item)
    for index,hole in enumerate(cfg.hole):
        if hole.mass<=0:continue
        item=dict(case='exact isolated seed',hole_index=index);output['controls'].append(item)
        sequence(lambda x,h=hole:backend.seed(h,x,cfg.conformal_choice),item)
    for record in records:
        config,unknowns=restore(backend,record,a.raw_directory)
        verify_saved_plan(config,record,a.raw_directory)
        new_xyz,new_roles,new_steps=sampling_plan(config,record,a.region,factors)
        if not (np.array_equal(xyz,new_xyz) and np.array_equal(roles,new_roles)
                and all(np.array_equal(s,t) for s,t in zip(steps,new_steps))):
            raise ValueError('every calibration grid must use identical points and steps')
        item=dict(resolution=record['resolution'],source_stopping_verified=record.get('stopping_verified',False))
        output['binary_grids'].append(item)
        with backend.create_sampler(config) as sampler:
            sampler.set_unknowns(unknowns);sequence(sampler.sample,item)
    output['completed']=True;save();print(target,flush=True)
    return 0


if __name__=='__main__':raise SystemExit(main())
