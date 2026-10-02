"""Serial matched fresh solves in separately verified native-library processes."""
import argparse,json,os,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
from hispid import Backend,Config,Hole
from configs import as_dict
from checkpoints import ROOT
from run_validation import points


def worker(library,path):
    b=Backend(library)
    evidence=json.loads((ROOT/'validation/polar_sequence_covariance.json').read_text())
    source=evidence['records'][0];c=b.config()
    for name,_ in Config._fields_:
        value=source['config'][name]
        if name=='hole':
            for h,hole in enumerate(value):c.hole[h]=Hole(**hole)
        elif hasattr(getattr(c,name),'_length_'):getattr(c,name)[:]=value
        else:setattr(c,name,value)
    xyz,near,bulk=points(c);xyz=xyz[:near+bulk]
    start=time.monotonic()
    with b.create(c) as s:
        diagnostic=s.solve();unknowns=s.unknowns();sample=s.sample(xyz)
        residual=s.residual(unknowns)
    report=dict(library_sha256=b.library_sha256(),loaded_image_verified=True,config=as_dict(c),
                diagnostics=diagnostic,whole_seconds=time.monotonic()-start,
                max_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
                initial_guess='zero',cpu_threads=1)
    np.savez_compressed(path,unknowns=unknowns,residual=residual,**sample)
    path.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report,flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--old-library');p.add_argument('--new-library')
    p.add_argument('--output',default='validation/preconditioner_reuse_benchmark.json')
    p.add_argument('--worker-library');p.add_argument('--worker-output');a=p.parse_args()
    if a.worker_library:worker(a.worker_library,Path(a.worker_output));return 0
    if not a.old_library or not a.new_library:p.error('both libraries required')
    out=ROOT/a.output;raw=ROOT/'validation/raw'/out.stem;raw.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    paths=[]
    for side,library in [('old',a.old_library),('new',a.new_library)]:
        path=raw/(side+'.npz');paths.append(path)
        subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker-library',str(Path(library).resolve(strict=True)),
                        '--worker-output',str(path)],cwd=ROOT,env=env,check=True,timeout=1200)
    old,new=[json.loads(path.with_suffix('.json').read_text()) for path in paths]
    with np.load(paths[0]) as x,np.load(paths[1]) as y:
        if set(x.files)!=set(y.files):raise ValueError('array inventories differ')
        differences={key:float(np.max(abs(x[key]-y[key]))) for key in x.files}
    counts_equal=all(old['diagnostics'][key]==new['diagnostics'][key] for key in ('status','newton_iterations','krylov_iterations'))
    result=dict(old=old,new=new,array_absolute_differences=differences,iteration_counts_identical=counts_equal,
                solve_speedup=old['diagnostics']['seconds']/new['diagnostics']['seconds'],
                isolated_processes=True,loaded_images_verified=True,
                passed=bool(old['config']==new['config'] and counts_equal and old['diagnostics']['status']==0 and max(differences.values())==0))
    result['note']='Same full free data/grid/tolerance/restart, zero initial guesses, one thread. Compare final iterates, physical fields and residual bitwise, and Newton/Krylov counts. Original binary acceptance records remain unchanged.'
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
    return 0 if result['passed'] else 1

if __name__=='__main__':raise SystemExit(main())
