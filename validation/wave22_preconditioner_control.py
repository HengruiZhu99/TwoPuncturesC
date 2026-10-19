"""Bounded algebraic control; the coarse grid does not qualify wave resolution."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Backend,Wave22
from kerr_wave_case import config


def worker(library,output):
    b=Backend(str(library.resolve()));c=config(b,32,.9999)
    c.n[:]=(32,64,16);c.hole[1].center[:]=(0,0,-16)
    c.krylov_restart=80;c.max_krylov=600;c.max_newton=16
    w=Wave22(amplitude=.001,omega=.9760972454792889,cycles=4,r0=12,cutoff_min=2,cutoff_max=4)
    rng=np.random.default_rng(473);u=rng.normal(0,1e-9,4*np.prod(c.n));d=rng.normal(0,1e-9,len(u))
    dirs=np.array([[.3,.4,.5],[-.7,.2,.3],[.2,-.8,.4]])
    dirs/=np.linalg.norm(dirs,axis=1)[:,None]
    xyz=(np.array([.12,.15,.19,.28,.35,1.,12.,24.,38.])[:,None,None]*dirs[None,:,:]).reshape(-1,3)
    with b.create(c,wave=w) as s:
        residual=s.residual(u);jvp=s.jvp(u,d)
        diag=s.solve(krylov='gmres',linear_rtol=.1)
        fields=s.sample(xyz)
        np.savez_compressed(output.with_suffix('.npz'),residual=residual,jvp=jvp,xyz=xyz,
                            gamma=fields['gamma'],Kij=fields['Kij'],unknowns=s.unknowns())
        result=dict(diagnostics=diag,linear_history=s.linear_history(),resolved=s.resolved_options,
                    library_sha256=b.library_sha256(),grid=list(c.n),nonlinear_tolerance=c.tolerance)
    output.write_text(json.dumps(result,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--old-library',type=Path);p.add_argument('--new-library',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--worker-library',type=Path,help=argparse.SUPPRESS)
    a=p.parse_args()
    if a.worker_library:
        worker(a.worker_library,a.output);return
    if not a.old_library or not a.new_library:p.error('both images required')
    a.output.mkdir(parents=True,exist_ok=False)
    for name,lib in [('old',a.old_library),('new',a.new_library)]:
        with (a.output/f'{name}.log').open('w') as log:
            subprocess.run([sys.executable,__file__,'--worker-library',str(lib.resolve()),'--output',str(a.output/f'{name}.json')],check=True,stdout=log,stderr=log)
    x=np.load(a.output/'old.npz');y=np.load(a.output/'new.npz')
    old=json.loads((a.output/'old.json').read_text());new=json.loads((a.output/'new.json').read_text())
    identical={k:bool(np.array_equal(x[k],y[k])) for k in ('residual','jvp')}
    result=dict(purpose=__doc__,operator_bitwise_identical=identical,
                all_native_converged=bool(old['diagnostics']['converged'] and new['diagnostics']['converged']),
                final_sample_max_difference={k:float(np.max(abs(x[k]-y[k]))) for k in ('gamma','Kij')},
                krylov_iterations={k:v['diagnostics']['krylov_iterations'] for k,v in [('old',old),('new',new)]},
                solve_seconds={k:v['diagnostics']['seconds'] for k,v in [('old',old),('new',new)]},
                driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (a.output/'comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
    if not all(identical.values()) or not result['all_native_converged']:raise SystemExit(1)


if __name__=='__main__':main()
