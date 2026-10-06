"""One small CPU/device control, not a repeated physical/backend matrix."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
from hispid import Backend,Hole
from execution import name,statistics,device_description

def run(path,geometries=('host','execution'),check_seeds=True):
    b=Backend(str(path.resolve()));c=b.config();c.n[:]=[8,12,4];c.seed_family='trumpet_r0_m'
    c.hole[0]=Hole(.5,center=(3,0,0),spin=(.02,-.03,.1),velocity=(.02,.08,-.03))
    c.hole[1]=Hole(.5,center=(-3,0,0),spin=(-.03,.02,.08),velocity=(-.02,-.08,.03))
    c.far_radius=0;c.inner_flatten=0;c.inner_min[:]=[0,0];c.inner_max[:]=[0,0]
    x=np.array([[1.2,.7,-.4],[2,-1,.8],[4,.3,2.]])
    seed=[]
    for h in ([Hole(1,spin=(0,0,.99)),Hole(1,velocity=(np.sqrt(.99),0,0)),c.hole[0]] if check_seeds else []):
        ref=b.seed(h,x,seed_family=c.seed_family)
        dev=b.seed(h,x,execution='kokkos',seed_family=c.seed_family)
        error=max(float(np.max(abs(dev[k]-ref[k])/(1+abs(ref[k])))) for k in ref)
        seed.append(error)
    with b.create(c) as reference:
        rng=np.random.default_rng(732);u=rng.normal(0,1e-5,reference.size);v=rng.normal(0,1e-3,reference.size)
        F=reference.residual(u);J=reference.jvp(u,v)
    rows=[]
    for geometry in geometries:
        start=time.monotonic()
        with b.create(c,execution='kokkos',geometry=geometry) as s:
            f=s.residual(u);j=s.jvp(u,v)
            rows.append(dict(geometry=geometry,residual_relative_l2=float(np.linalg.norm(f-F)/np.linalg.norm(F)),
                jvp_relative_l2=float(np.linalg.norm(j-J)/np.linalg.norm(J)),seconds=time.monotonic()-start))
    return dict(kind='trumpet_small_device_control',backend=name(b.lib),device=device_description(b.lib),
                library_sha256=b.library_sha256(),collocation_maps=b.parameterization_maps(),
                seed_control_performed=check_seeds,driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seed_scaled_errors=seed,rows=rows,statistics=statistics(b.lib),
                passed=max(seed,default=0)<1e-10 and all(max(r['residual_relative_l2'],r['jvp_relative_l2'])<1e-7 for r in rows),
                binary_acceptance=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--library',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--geometry',nargs='+',choices=['host','execution'],default=['host','execution']);ap.add_argument('--skip-seeds',action='store_true');a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=run(a.library,a.geometry,not a.skip_seeds);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
