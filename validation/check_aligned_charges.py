"""Independent Cartesian-FD fluxes on the native aligned sphere grid.

This is an integration/API control only. It never solves or accepts a binary.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from hispid import Backend
from checkpoints import restore_equivalent,select_record
from physical import charges

def main():
    p=argparse.ArgumentParser();p.add_argument('--library',required=True)
    p.add_argument('--proof',required=True);p.add_argument('--output',required=True)
    p.add_argument('--case',default='moderate_default_polar2_r128')
    p.add_argument('--radii',default='40,200,1000');p.add_argument('--quadrature',default='16:32')
    a=p.parse_args();backend=Backend(a.library);record=select_record(a.case)
    cfg,unknowns=restore_equivalent(backend,record,json.loads(Path(a.proof).read_text()))
    direction=np.array(cfg.hole[0].center)-np.array(cfg.hole[1].center)
    direction/=np.linalg.norm(direction);reference=np.eye(3)[np.argmin(abs(direction))]
    tangent=reference-direction*(reference@direction);tangent/=np.linalg.norm(tangent)
    frame=np.column_stack((direction,tangent,np.cross(direction,tangent)))
    axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis)
    W=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    Q=np.eye(3)+np.sin(.73)*W+(1-np.cos(.73))*(W@W);offset=np.array([.7,-.2,.4])
    nt,np_=map(int,a.quadrature.split(':'));rows=[]
    # FD truncation and summation use separate implementations. Bound declared
    # before this control; no existing binary/seed acceptance limits are altered.
    result=dict(case=a.case,resolution=record['resolution'],library_sha256=backend.library_sha256(),
        checkpoint_source_library_sha256=record['library_sha256'],quadrature=[nt,np_],
        polar_frame=frame.tolist(),charge_absolute_tolerance=1e-8,
        centered_covariance_absolute_tolerance=1e-8,records=rows,passed=False,binary_acceptance=False)
    with backend.create_sampler(cfg) as solution:
        solution.set_unknowns(unknowns)
        def transformed(x):
            fields=solution.sample((np.asarray(x)-offset)@Q)
            for name in ('gamma','Kij'):
                fields[name]=np.einsum('ik,nkl,jl->nij',Q,fields[name].reshape(-1,3,3),Q).reshape(-1,9)
            return fields
        for radius in map(float,a.radii.split(',')):
            start=time.monotonic();native=solution.charges(radius,ntheta=nt,nphi=np_)
            independent=charges(solution.sample,radius,ntheta=nt,nphi=np_,polar_frame=frame)
            rotated=charges(transformed,radius,center=offset,ntheta=nt,nphi=np_,polar_frame=Q@frame)
            expected=np.r_[independent[0],Q@independent[1:4],Q@independent[4:7]]
            difference=float(np.max(abs(native-independent)));covariance=float(np.max(abs(rotated-expected)))
            row=dict(radius=radius,native_EPJ=native.tolist(),independent_FD_EPJ=independent.tolist(),
                charge_difference_linf=difference,centered_covariance_difference_linf=covariance,
                seconds=time.monotonic()-start,passed=bool(difference<1e-8 and covariance<1e-8))
            rows.append(row);Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(row,flush=True)
    result['passed']=all(r['passed'] for r in rows);Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    if not result['passed']:raise SystemExit(1)

if __name__=='__main__':main()
