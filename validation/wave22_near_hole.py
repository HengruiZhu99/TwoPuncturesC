"""Fixed off-grid correction/geometry differences across the near-hole transition.

This is a solution-difference diagnostic, not a vacuum-constraint test inside
an intentionally modified interior. No residual or stopping norm is rescaled.
"""
import argparse
import copy
import ctypes as C
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from checkpoint_export import read_checkpoint
from hispid import Backend
from parallel_observer import ParallelObserver
from physical import constraints,norms


def signature(c):
    out={}
    for name,_ in c._fields_:
        if name in ('n','max_newton','max_krylov','krylov_restart','memory_limit_mib'):continue
        v=getattr(c,name)
        out[name]=[[h.mass,*h.center,*h.spin,*h.velocity] for h in v] if name=='hole' else (list(v) if isinstance(v,C.Array) else v)
    out.update(seed_family=c.seed_family,wave22=c.wave22.as_dict())
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library',type=Path,required=True)
    p.add_argument('--checkpoints',type=Path,nargs=3,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--threads',type=int,default=16)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    b=Backend(str(a.library.resolve()));reference=None;rows=[];data=[];archive={};seed_floor=None
    radius=np.unique(np.r_[np.linspace(.08,.4,97),.1,.15,.17451721,.2,.25353544551647955,.2939671,.30072])
    mu=np.array([-.95,-.6,0.,.6,.95]);phi=np.pi/8
    rays=np.stack((np.sqrt(1-mu*mu)*np.cos(phi),np.sqrt(1-mu*mu)*np.sin(phi),mu),axis=-1)
    xyz=(radius[:,None,None]*rays[None,:,:]).reshape(-1,3)
    exterior_radius=np.array([.21,.24,1.05*.25353544551647955,.2939671,.35,.4])
    exterior_xyz=(exterior_radius[:,None,None]*rays[None,:,:]).reshape(-1,3)
    exterior_steps=.006*np.linalg.norm(exterior_xyz,axis=1)
    # Even the innermost FD6 stencil stays outside the modified r<=.2 ball.
    if np.min(np.linalg.norm(exterior_xyz,axis=1)-3*exterior_steps)<=.2:raise ValueError('constraint stencil intersects modified interior')
    regions=dict(flat_interior=np.repeat(radius<.1,len(mu)),interior_transition=np.repeat((radius>=.1)&(radius<=.2),len(mu)),near_horizon=np.repeat(radius>.2,len(mu)))
    for index,path in enumerate(a.checkpoints):
        c,u,meta=read_checkpoint(path);case=signature(c)
        if rows and c.n[0]<=rows[-1]['n'][0]:raise ValueError('increasing radial resolutions required')
        if meta['acceptance']=='diagnostic':raise ValueError('strict native convergence required')
        if meta['parameterization']!=b.parameterization():raise ValueError('basis mismatch')
        if reference is not None and case!=reference:raise ValueError('formulation or tolerance changed')
        reference=case
        if (c.seed_family!='liu' or c.hole[0].mass!=1 or any(c.hole[0].center) or c.hole[1].mass!=0
                or list(c.hole[0].spin)!=[0.,0.,.9999] or any(c.hole[0].velocity)
                or not c.inner_flatten or c.inner_min[0]!=.1 or c.inner_max[0]!=.2):
            raise ValueError('this fixed near-hole observer requires the centered, unboosted mass-one Liu chi=.9999 case with transition [.1,.2]')
        with b.create_sampler(c) as s:
            s.set_unknowns(u)
            with ParallelObserver(s,a.threads) as observer:
                q=observer(xyz);exterior=norms(constraints(observer,exterior_xyz,exterior_steps,6))
            if seed_floor is None:
                unperturbed=copy.copy(c);unperturbed.n[:]=(8,8,8);del unperturbed.wave22
                with b.create_sampler(unperturbed) as seed,ParallelObserver(seed,a.threads) as observer:
                    seed_floor=norms(constraints(observer,exterior_xyz,exterior_steps,6))
        row=dict(n=list(c.n),checkpoint=meta,exterior_constraints=exterior)
        if data:
            row['difference_from_previous']={region:{key:dict(linf=float(np.max(abs(q[key][mask]-data[-1][key][mask]))),rms=float(np.sqrt(np.mean((q[key][mask]-data[-1][key][mask])**2)))) for key in ('correction','gamma','Kij')} for region,mask in regions.items()}
        archive.update({f'level{index}_{key}':q[key] for key in ('correction','gamma','Kij')})
        rows.append(row);data.append(q)
    decreasing={region:{key:rows[2]['difference_from_previous'][region][key]['linf']<rows[1]['difference_from_previous'][region][key]['linf'] for key in ('correction','gamma','Kij')} for region in regions}
    np.savez_compressed(a.output/'profiles.npz',radius=radius,mu=mu,xyz=xyz,**archive)
    constraint_decreases={k:all(rows[i+1]['exterior_constraints'][k]<rows[i]['exterior_constraints'][k] for i in (0,1)) for k in ('H_rms','M_rms')}
    required_regions=('interior_transition','near_horizon')
    passed=all(constraint_decreases.values()) and all(decreasing[r][k] for r in required_regions for k in ('correction','gamma','Kij'))
    result=dict(purpose=__doc__,library_sha256=b.library_sha256(),driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),case=reference,rows=rows,
                exterior_constraint_observer=dict(xyz=exterior_xyz.tolist(),steps=exterior_steps.tolist(),order=6,norm='unweighted fixed-point RMS',seed_calibration=seed_floor,note='Exterior only, includes equator; raw calibration is retained without subtraction.'),
                exterior_constraint_rms_decreases=constraint_decreases,required_solution_regions=list(required_regions),passed=bool(passed),
                points=len(xyz),radius_range=[float(radius[0]),float(radius[-1])],phi=float(phi),mu=mu.tolist(),linf_differences_decrease=decreasing,
                note='Raw absolute solution differences. Modified-interior physical constraints are not required to vanish; no continuum error bound is inferred.',profiles_sha256=hashlib.sha256((a.output/'profiles.npz').read_bytes()).hexdigest())
    (a.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
    if not passed:raise SystemExit('near-hole Cauchy differences or exterior constraint RMS failed to decrease')


if __name__=='__main__':main()
