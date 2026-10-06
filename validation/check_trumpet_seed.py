"""Focused independent physical/derivative check of the development seed probe.

Uses raw physical tensors and the existing independent Cartesian observer.
No performance measurement or multi-backend matrix.
"""
import argparse
import ctypes as C
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
from hispid import Hole, PTR, ptr
from physical import constraints, norms

def run(library, step_scale=1.):
    lib=C.CDLL(str(library.resolve()))
    lib.trumpet_probe.argtypes=[C.POINTER(Hole),C.c_int,PTR,PTR,PTR,PTR,PTR]
    lib.trumpet_probe.restype=C.c_int
    rows=[]
    for label,spin,velocity in [('static',(0,0,0),(0,0,0)),
                                ('generic',(.2,-.3,.7),(.25,.1,-.2)),
                                ('gamma10',(0,0,0),(np.sqrt(.99),0,0))]:
        if not np.isfinite(step_scale) or step_scale<=0: raise ValueError("positive finite step scale required")
        hole=Hole(1,spin=spin,velocity=velocity)
        def sample(x):
            x=np.ascontiguousarray(x,dtype=float).reshape(-1,3);n=len(x)
            g=np.empty((n,3,3));k=np.empty_like(g)
            dg=np.empty((n,3,3,3));dk=np.empty_like(dg)
            status=lib.trumpet_probe(C.byref(hole),n,ptr(x),ptr(g),ptr(k),ptr(dg),ptr(dk))
            if status:raise RuntimeError(f'{label}: rejected seed status {status}')
            return dict(gamma=g,Kij=k,dgamma=dg,dKij=dk,attenuation=np.ones(n))
        x=np.array([[1.2,.7,-.4],[2,-1,.8],[4,.3,2.]])
        if label=='gamma10':x[:,0]/=10
        levels=[]
        for h in [.008,.004,.002]:
            h*=step_scale
            if label=='gamma10':h/=10
            levels.append(dict(step=h,**norms(constraints(sample,x,h))))
        center=sample(x);derivative_error={}
        for field,key in [('gamma','dgamma'),('Kij','dKij')]:
            errs=[]
            for h in [.004,.002,.001]:
                h*=step_scale
                if label=='gamma10':h/=10
                approx=np.empty_like(center[key])
                for d in range(3):
                    shift=np.zeros(3);shift[d]=h
                    approx[:,d]=(sample(x-2*shift)[field]-8*sample(x-shift)[field]
                        +8*sample(x+shift)[field]-sample(x+2*shift)[field])/(12*h)
                errs.append(float(np.max(np.abs(approx-center[key])/(1+np.abs(center[key])))))
            derivative_error[field]=errs
        rows.append(dict(case=label,constraints=levels,derivative_scaled_linf=derivative_error,
                         fine_constraint_passed=levels[-1]['H_rms']<1e-7 and levels[-1]['M_rms']<1e-7))
    root=Path(__file__).resolve().parents[1]
    return dict(kind='focused_development_seed_check',rows=rows,
                full_seed_acceptance=False,binary_acceptance=False,
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                    [library.resolve(),Path(__file__).resolve(),root/'validation/physical.py',
                     root/'tests/trumpet_probe.cpp',root/'src/HiSpID_trumpet.hpp']})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--library',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--step-scale',type=float,default=1.);args=ap.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    result=run(args.library,args.step_scale);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['rows'],indent=2))
