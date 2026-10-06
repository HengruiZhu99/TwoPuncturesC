"""One independent stencil sensitivity check plus missing coarse ADM charges."""
import argparse,json,sys,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'python'))
from hispid import Backend
from checkpoint_export import read_checkpoint
from physical import constraints,norms,extrapolate

def run(library,pilot):
    b=Backend(str(library.resolve()));c,u,meta=read_checkpoint(pilot/'diagnostic.checkpoint')
    if meta['source_library_sha256']!=b.loaded_sha256:raise ValueError('source image mismatch')
    with np.load(pilot/'physical.npz') as raw:
        x,step,H,M=raw['xyz'],raw['step'],raw['H'],raw['Mnorm']
    old=json.loads((pilot/'result.json').read_text());n=old['physical']['near']['count'];nb=old['physical']['bulk']['count']
    with b.create_sampler(c) as s:
        s.set_unknowns(u);r=constraints(s.sample,x[:n+nb],step[:n+nb]/2)
        radii=[256.,512.,1024.];q=[s.charges(R,ntheta=2*c.n[1],nphi=32).tolist() for R in radii]
    return dict(checkpoint=meta,driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        half_step_physical=dict(near=norms(r,np.arange(n+nb)<n),bulk=norms(r,np.arange(n+nb)>=n)),
        H_stencil_change_linf=float(np.max(abs(r['H']-H[:n+nb]))),Mnorm_stencil_change_linf=float(np.max(abs(r['Mnorm']-M[:n+nb]))),
        charge_radii=radii,charges=q,charges_extrapolated=extrapolate(radii,q).tolist(),binary_acceptance=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--library',type=Path,required=True);ap.add_argument('--pilot',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=run(a.library,a.pilot);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
