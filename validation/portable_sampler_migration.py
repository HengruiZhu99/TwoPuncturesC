"""Prove a particular portable checkpoint's CPU sampler migration.

The coordinator loads no native image. Producer and consumer replay precisely
the same coefficients in separate processes, including physical derivatives
on trial horizon points. This proves sampling compatibility only: it neither
re-solves the PDE nor transfers a checkpoint's physical acceptance.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from checkpoint_export import read_checkpoint
from hispid import Backend,Point
from run_validation import points
from execution import name as execution_name
from native_loader import loaded_kokkos_images

TOLERANCE=1e-12
FIELDS=tuple(name for name,_ in Point._fields_)+('dgamma',)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def witness_points(config):
    base,_,_=points(config,True)
    mu,_=np.polynomial.legendre.leggauss(12)
    phi=2*np.pi*(np.arange(24)+.37)/24
    mu,phi=np.meshgrid(mu,phi,indexing='ij')
    directions=np.c_[np.sqrt(1-mu.ravel()**2)*np.cos(phi.ravel()),
                     np.sqrt(1-mu.ravel()**2)*np.sin(phi.ravel()),mu.ravel()]
    directions=np.r_[directions,np.eye(3),-np.eye(3)]
    trials=[]
    for hole in config.hole:
        if hole.mass<=0:continue
        spin=np.asarray(hole.spin);v=np.asarray(hole.velocity)
        chi=np.linalg.norm(spin)/hole.mass**2
        if not 0<=chi<1 or not v@v<1:raise ValueError('invalid subextremal seed')
        radius=.5*hole.mass*np.sqrt(1-chi**2)/np.sqrt(1+(directions@v)**2/(1-v@v))
        trials.extend(np.asarray(hole.center)+factor*radius[:,None]*directions for factor in (.8,1,1.2))
    return np.ascontiguousarray(np.concatenate([base,*trials,100*directions,1000*directions]))


def worker(args):
    cfg,unknowns,checkpoint=read_checkpoint(args.checkpoint)
    backend=Backend(str(Path(args.worker_library).resolve(strict=True)))
    if backend.library_sha256()!=args.worker_sha:raise ValueError('worker library mismatch')
    if backend.parameterization()!=checkpoint['parameterization']:
        raise ValueError('different continuous basis cannot be a sampler migration')
    xyz=witness_points(cfg)
    with backend.create_sampler(cfg) as sampler:
        sampler.set_unknowns(unknowns)
        values=sampler.sample_with_derivatives(xyz)
    if digest(args.checkpoint)!=checkpoint['file_sha256']:raise ValueError('checkpoint changed during sampler witness')
    if set(values)!=set(FIELDS):raise ValueError('unexpected field inventory')
    if not all(np.isfinite(v).all() for v in values.values()):raise ValueError('nonfinite sampler witness')
    output=Path(args.worker_output)
    if output.exists() or output.with_suffix('.json').exists():raise FileExistsError(output)
    np.savez_compressed(output,xyz=xyz,**values)
    output.with_suffix('.json').write_text(json.dumps(dict(
        library_path=str(backend.path),library_sha256=backend.library_sha256(),
        dependency_images=backend.dependency_images,checkpoint=checkpoint,
        compiled_execution=execution_name(backend.lib),runtime_images=loaded_kokkos_images(),
        parameterization=backend.parameterization(),maps=backend.parameterization_maps(),
        point_count=len(xyz),artifact_sha256=digest(output)),indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint',required=True)
    parser.add_argument('--producer-library');parser.add_argument('--consumer-library')
    parser.add_argument('--output')
    parser.add_argument('--timeout',type=int,default=3600)
    parser.add_argument('--worker-library');parser.add_argument('--worker-sha');parser.add_argument('--worker-output')
    args=parser.parse_args()
    if args.worker_library:
        if not args.worker_sha or not args.worker_output:parser.error('worker SHA and output required')
        worker(args);return 0
    if not args.producer_library or not args.consumer_library or not args.output:
        parser.error('producer, consumer and output required')
    if args.timeout<1:parser.error('positive worker timeout required')
    checkpoint_path=Path(args.checkpoint).resolve(strict=True)
    _,_,checkpoint=read_checkpoint(checkpoint_path)
    libraries=[Path(p).resolve(strict=True) for p in (args.producer_library,args.consumer_library)]
    hashes=[digest(p) for p in libraries]
    if hashes[0]!=checkpoint['source_library_sha256'] or hashes[0]==hashes[1]:
        raise ValueError('bound producer and distinct consumer required')
    output=Path(args.output).resolve();output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():raise FileExistsError(output)
    raw=output.parent/'raw'/output.stem;raw.mkdir(parents=True,exist_ok=False)
    result=dict(schema='hispid_portable_sampler_migration_v1',purpose='sampling_only_no_acceptance_transfer',
        checkpoint=checkpoint,producer_library_sha256=hashes[0],consumer_library_sha256=hashes[1],
        isolated_processes=True,loaded_images_verified=False,
        criteria=dict(field_scaled_linf=TOLERANCE),workers=[],differences={},passed=False)
    def save():output.write_text(json.dumps(result,indent=2)+'\n')
    save()
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    for side,library,sha in zip(('producer','consumer'),libraries,hashes):
        if digest(checkpoint_path)!=checkpoint['file_sha256']:raise ValueError('checkpoint changed before worker')
        if any(digest(p)!=s for p,s in zip(libraries,hashes)):raise ValueError('bound library changed before worker')
        artifact=raw/(side+'.npz')
        command=[sys.executable,str(Path(__file__).resolve()),'--checkpoint',str(checkpoint_path),
                 '--worker-library',str(library),'--worker-sha',sha,'--worker-output',str(artifact)]
        record=dict(side=side,command=command,passed=False)
        result['workers'].append(record);save()
        with (raw/(side+'.log')).open('w') as log:
            try:run=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout)
            except subprocess.TimeoutExpired:record['returncode']='timeout';save();return 1
        record['returncode']=run.returncode
        if run.returncode:save();return 1
        record.update(json.loads(artifact.with_suffix('.json').read_text()),artifact=str(artifact))
        if digest(checkpoint_path)!=checkpoint['file_sha256']:raise ValueError('checkpoint changed during worker')
        record['passed']=bool(record['library_sha256']==sha and record['checkpoint']==checkpoint
                              and record['artifact_sha256']==digest(artifact))
        save()
    old,new=result['workers']
    result['pure_reference_consumer']=bool(new['compiled_execution']=='reference' and not new['runtime_images']
        and not any('kokkos' in Path(p).name.lower() for p in new['dependency_images']))
    result['loaded_images_verified']=all(w['passed'] for w in result['workers'])
    result['parameterization']=old['parameterization'];result['maps']=old['maps']
    with np.load(old['artifact']) as a,np.load(new['artifact']) as b:
        if set(a.files)!=set(FIELDS)|{'xyz'} or set(a.files)!=set(b.files):raise ValueError('witness inventory differs')
        if not np.array_equal(a['xyz'],b['xyz']):raise ValueError('witness points differ')
        result['point_count']=len(a['xyz'])
        for field in FIELDS:
            x,y=a[field],b[field]
            if x.shape!=y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('invalid field witness')
            result['differences'][field]=float(np.max(np.abs(x-y)/(1+np.abs(x))))
        result['positive_metric']=bool(np.linalg.eigvalsh(a['gamma'].reshape(-1,3,3)).min()>0
                                      and np.linalg.eigvalsh(b['gamma'].reshape(-1,3,3)).min()>0)
    result['passed']=bool(result['loaded_images_verified'] and result['positive_metric'] and result['pure_reference_consumer']
        and old['parameterization']==new['parameterization'] and old['maps']==new['maps']
        and max(result['differences'].values())<=TOLERANCE)
    result['note']='Identical coefficients and exact basis/maps; separate loaded-image witnesses, off-grid and trial-horizon fields/derivatives. This empirical sampler comparison does not establish physical constraints, exact PDE equivalence, or horizon accuracy. Source acceptance is preserved.'
    if digest(checkpoint_path)!=checkpoint['file_sha256']:raise ValueError('checkpoint changed before qualification')
    for worker in result['workers']:
        if (digest(worker['library_path'])!=worker['library_sha256']
            or digest(worker['artifact'])!=worker['artifact_sha256']
            or any(digest(p)!=s for p,s in (worker['dependency_images']|worker['runtime_images']).items())):
            raise ValueError('sampler image or raw witness changed before qualification')
    save();print(json.dumps(result['differences']),flush=True)
    return 0 if result['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
