"""Replay a retained iterate on a pure CPU reference, without solving again."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from checkpoint_export import read_checkpoint
from hispid import Backend
from configs import as_dict
from physical import constraints, norms
from native_loader import loaded_kokkos_images
from execution import name as execution_name


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('checkpoint','row-receipt','reference-library','output'):
        p.add_argument('--'+key,required=True)
    a=p.parse_args();root=Path(a.output).resolve();root.mkdir(parents=True,exist_ok=False)
    checkpoint=Path(a.checkpoint).resolve(strict=True)
    receipt=Path(a.row_receipt).resolve(strict=True)
    row=json.loads(receipt.read_bytes())['record']
    cfg,values,cp=read_checkpoint(checkpoint)
    if cp['source_library_sha256']!=row['library_sha256'] or as_dict(cfg)!=row['config']:
        raise ValueError('checkpoint and retained solve configuration differ')
    raw=row['raw_artifact_sha256']
    fields=next(Path(x) for x in raw if x.endswith('_'+str(cfg.n[2])+'.npz'))
    nodes=next(Path(x) for x in raw if x.endswith('_collocation.npz'))
    bound={str(checkpoint):digest(checkpoint),str(receipt):digest(receipt),**raw}
    bound[str(Path(__file__).resolve())]=digest(__file__)
    for module in ('checkpoint_export','hispid','configs','physical','native_loader','execution'):
        file=Path(__import__(module).__file__).resolve();bound[str(file)]=digest(file)
    def verify():
        if any(digest(file)!=sha for file,sha in bound.items()):
            raise ValueError('bound replay source, image or artifact changed')
    verify()
    for file in (fields,next(Path(x) for x in raw if x.endswith('_solve.npz'))):
        with np.load(file) as saved:
            reference_values=saved['unknowns']
            if (values.shape!=reference_values.shape or values.dtype!=reference_values.dtype
                or not np.array_equal(values.view(np.uint64),reference_values.view(np.uint64))):
                raise ValueError('checkpoint and SHA-bound retained iterate differ')
    result=dict(purpose='fixed_iterate_CPU_equation_and_reconstruction_diagnostic',
        completed=False,solve_performed=False,physical_acceptance=False,
        bound_artifacts_sha256=bound,FD_records=[])
    def save():
        verify();(root/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    save()
    backend=Backend(a.reference_library)
    if (loaded_kokkos_images() or execution_name(backend.lib)!='reference'
        or backend.parameterization()!=cp['parameterization']
        or backend.parameterization_maps()!=row['collocation_maps']
        or backend.residual_scaling()!=row['residual_scaling']):
        raise ValueError('pure reference and identical continuous basis, maps and norm required')
    result['checkpoint_matches_bound_iterate']=True
    result['reference_basis_maps_norm_verified']=True
    bound.update({str(backend.path):backend.library_sha256(),**backend.dependency_images})
    result['reference_images']={str(backend.path):backend.library_sha256(),**backend.dependency_images}
    result['source_diagnostics']=row['diagnostics'];save()
    with np.load(fields) as f:
        nnear=row['near_sample_count'];nbulk=row['bulk_sample_count']
        worst=nnear+np.argsort(np.abs(f['H'][nnear:nnear+nbulk]))[-3:]
        worst_xyz=f['points'][worst].copy()
    with np.load(nodes) as f:
        source={key:f[key].copy() for key in f.files}
    with backend.create(cfg,execution='reference',geometry='host') as s:
        s.set_unknowns(values)
        residual=s.residual(values).reshape(-1,4)
        current=s.equation_samples()
        result['CPU_weighted_linf']=np.max(np.abs(residual),axis=0).tolist()
        result['CPU_internal_stopping_verified']=bool(np.isfinite(residual).all() and np.max(np.abs(residual))<=cfg.tolerance)
        result['collocation_differences']={key:float(np.max(np.abs(current[key]-source[key])/(1+np.abs(source[key])))) for key in source}
        eligible=np.flatnonzero(current['attenuation']==1)
        chosen=np.unique([eligible[np.argmin(np.sum((current['xyz'][eligible]-x)**2,axis=1))] for x in worst_xyz])
        node_xyz=current['xyz'][chosen]
        sampled=s.sample(node_xyz)
        result['node_indices']=chosen.tolist();result['node_xyz']=node_xyz.tolist()
        result['native_node_psi']=current['psi'][chosen].tolist()
        result['sampled_node_psi']=sampled['psi'].tolist()
        result['node_psi_scaled_difference']=float(np.max(np.abs(sampled['psi']-current['psi'][chosen])/(1+np.abs(current['psi'][chosen]))))
        result['worst_offgrid_xyz']=worst_xyz.tolist();save()
        np.savez_compressed(root/'CPU_collocation.npz',weighted=residual,**current)
        query=np.r_[node_xyz,worst_xyz]
        steps=np.minimum(.002,.001*np.min([np.linalg.norm(query-np.asarray(h.center),axis=1) for h in cfg.hole if h.mass>0],axis=0))
        for factor in (2,1,.5):
            physical=constraints(s.sample,query,factor*steps)
            raw_path=root/('physical_step_'+str(factor)+'.npz')
            sampled_fields={'sampled_'+key:value for key,value in s.sample_with_derivatives(query).items()}
            np.savez_compressed(raw_path,xyz=query,steps=factor*steps,**physical,**sampled_fields)
            result['FD_records'].append(dict(step_factor=factor,
                collocation_nodes=norms(physical,np.arange(len(query))<len(node_xyz)),
                offgrid=norms(physical,np.arange(len(query))>=len(node_xyz)),
                raw_artifact=str(raw_path),raw_sha256=digest(raw_path)))
            save()
    result['output_artifacts_sha256']={str(f):digest(f) for f in root.glob('*.npz')}
    result['completed']=True;save()
    print(json.dumps({key:result[key] for key in ('CPU_weighted_linf','collocation_differences','node_psi_scaled_difference','FD_records')}),flush=True)
    return 0


if __name__=='__main__':raise SystemExit(main())
