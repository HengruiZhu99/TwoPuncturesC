"""Calibrate isolated-vacuum source cancellation in a build's actual norm."""
import argparse,json,time,hashlib,os
from pathlib import Path
import numpy as np
from hispid import Backend,Hole
from configs import as_dict
from native_loader import loaded_kokkos_images

def isolated_controls(extreme=False):
    spin=np.array([.2,-.3,.4]);boost=np.array([-.7,.2,.3])
    z=np.zeros(3)
    if extreme:
        speed=float(np.sqrt(.99))
        spin*=.99/np.linalg.norm(spin);boost*=speed/np.linalg.norm(boost)
        return [('Schwarzschild',z,z),('spin99',np.array([0,0,.99]),z),
                ('spin99_generic',spin,z),('gamma10',z,np.array([speed,0,0])),
                ('gamma10_generic',z,boost)]
    spin*=.95/np.linalg.norm(spin);boost*=.885/np.linalg.norm(boost)
    return [('Schwarzschild',z,z),('spin95',spin,z),('boost885',z,boost),
            ('combined_generic',spin,boost)]

def main():
    p=argparse.ArgumentParser();p.add_argument('--library',required=True);p.add_argument('--output',required=True)
    p.add_argument('--resolutions',default='32:32:16,64:64:28,80:160:28')
    p.add_argument('--extreme',action='store_true',help='fresh exact chi.99/Gamma10 source-cancellation controls')
    p.add_argument('--execution',choices=('reference','kokkos'),default='reference')
    p.add_argument('--geometry',choices=('host','execution'),default='host')
    p.add_argument('--threads',type=int,default=1);p.add_argument('--memory-mib',type=int,default=8192)
    p.add_argument('--seed-mass',type=float,default=1.)
    p.add_argument('--coordinate-separation',type=float,default=12.,help='extreme controls use the same two chart centers as the binary')
    a=p.parse_args();path=Path(a.output)
    if a.geometry=='execution' and a.execution!='kokkos':raise ValueError('execution geometry requires Kokkos')
    if path.exists():raise FileExistsError('preserve prior source-floor controls')
    if not 1<=a.threads<=16 or not 1<=a.memory_mib<=65536:raise ValueError('invalid concurrency or memory budget')
    if not np.isfinite([a.seed_mass,a.coordinate_separation]).all() or min(a.seed_mass,a.coordinate_separation)<=0:
        raise ValueError('finite positive mass and chart separation required')
    os.environ.update(OMP_NUM_THREADS=str(a.threads),OPENBLAS_NUM_THREADS='1',OMP_PROC_BIND='close',OMP_PLACES='cores')
    backend=Backend(a.library);records=[]
    from execution import select,name,concurrency,device_description
    select(backend.lib,a.execution,a.threads)
    device=device_description(backend.lib)
    if a.execution=='kokkos' and name(backend.lib)=='Cuda' and (not device or device['visible_count']!=1 or device['visible_ordinal']!=0):
        raise ValueError('one allocated GPU required')
    if a.extreme and backend.residual_scaling()!='sin3_alpha_beta':raise ValueError('actual cubic weighted norm required')
    frozen={str(backend.path):backend.library_sha256()}|backend.dependency_images|loaded_kokkos_images()
    def save():
        if any(hashlib.sha256(Path(file).read_bytes()).hexdigest()!=sha for file,sha in frozen.items()):
            raise ValueError('bound source-floor image changed')
        for row in records:
            artifact=row.get('raw_artifact')
            if artifact and hashlib.sha256(Path(artifact['path']).read_bytes()).hexdigest()!=artifact['sha256']:
                raise ValueError('retained source-floor arrays changed')
        path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(output,indent=2)+'\n')
    output=dict(schema='hispid_source_floor_v3',library_sha256=backend.library_sha256(),residual_scaling=backend.residual_scaling(),
        weighted_far_source_bound=1e-14,far_radius_minimum=100.,exact_correction=0.,records=records,passed=False,
        extreme_controls=a.extreme,execution=a.execution,geometry=a.geometry,compiled_execution=name(backend.lib),
        seed_mass=a.seed_mass,coordinate_separation=a.coordinate_separation,
        execution_concurrency=concurrency(backend.lib),device=device,bound_images=frozen,
        unknown_parameterization_id=backend.parameterization(),collocation_maps=backend.parameterization_maps(),
        acceptance=False,note='Exact isolated vacuum seeds with all attenuation disabled; zero corrections. The far source floor is measured before any binary in a new norm. This does not accept a binary.')
    save()
    for shape in map(lambda x:list(map(int,x.split(':'))),a.resolutions.split(',')):
        if len(shape)!=3 or min(shape)<4:raise ValueError('three supported grid dimensions required')
        controls=[(label,S,v,active) for label,S,v in isolated_controls(a.extreme)
                  for active in ((0,1) if a.extreme else (0,))]
        for label,S,v,active in controls:
            cfg=backend.config();cfg.n[:]=shape;cfg.memory_limit_mib=a.memory_mib;cfg.krylov_restart=32
            if a.extreme:
                for index,sign in enumerate((1,-1)):cfg.hole[index]=Hole(0,(sign*a.coordinate_separation/2,0,0))
                sign=1 if active==0 else -1
                cfg.hole[active]=Hole(a.seed_mass,(sign*a.coordinate_separation/2,0,0),S*a.seed_mass**2,-sign*v)
            else:
                cfg.hole[0]=Hole(a.seed_mass,spin=S*a.seed_mass**2,velocity=v)
                cfg.hole[1]=Hole(0,center=(-a.coordinate_separation/2,0,0))
            cfg.conformal_choice=0;cfg.inner_flatten=0;cfg.inner_min[:]=cfg.inner_max[:]=[0,0]
            cfg.omega[:]=[0,0];cfg.far_radius=0;start=time.monotonic()
            row=dict(case=label,active_hole=active,resolution=shape,config=as_dict(cfg),completed=False,passed=False)
            records.append(row);save()
            with backend.create(cfg,execution=a.execution,geometry=a.geometry) as solution:
                row['setup_statistics']=solution.setup_statistics()
                if a.geometry=='execution' and (not row['setup_statistics'] or
                        row['setup_statistics']['geometry_execution']!=1 or row['setup_statistics']['scalar_digits']!=53):
                    raise ValueError('execution-built double geometry witness required')
                unknowns=solution.unknowns()
                if not np.all(unknowns==0):raise ValueError('exact isolated control requires zero unknowns')
                weighted=solution.residual(unknowns).reshape(-1,4)
                sampled=solution.equation_samples()
            radius=np.linalg.norm(sampled['xyz'],axis=1);mask=radius>=100
            if not np.any(mask):raise ValueError('no far nodes')
            maxima=np.max(abs(weighted[mask]),axis=0)
            row.update(far_node_count=int(mask.sum()),far_radius_range=[float(radius[mask].min()),float(radius[mask].max())],
                weighted_far_source_linf=maxima.tolist(),physical_equivalent_far_linf=np.max(abs(sampled['physical_equivalent'][mask]),axis=0).tolist(),
                passed=bool(np.all(np.isfinite(maxima)) and maxima.max()<1e-14),completed=True,seconds=time.monotonic()-start)
            artifact=path.parent/'raw'/path.stem/(label+'_hole'+str(active)+'_'+'_'.join(map(str,shape))+'.npz')
            artifact.parent.mkdir(parents=True,exist_ok=True)
            if artifact.exists():raise FileExistsError('preserve source-floor arrays')
            np.savez_compressed(artifact,xyz=sampled['xyz'],weighted=weighted,
                physical_equivalent=sampled['physical_equivalent'],far_mask=mask,unknowns=unknowns)
            row['raw_artifact']=dict(path=str(artifact.resolve()),sha256=hashlib.sha256(artifact.read_bytes()).hexdigest())
            save();print(label,'hole',active,shape,maxima,row['passed'],flush=True)
    output['passed']=bool(records and all(r['passed'] and r['completed'] for r in records));save()
    if not output['passed']:raise SystemExit(1)

if __name__=='__main__':main()
