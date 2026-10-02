"""Fixed exterior axis/near-axis Cartesian constraints on bound checkpoints.

This supplements ordinary off-grid directions; sampled checks do not certify
a continuous supremum over the whole exterior.
"""
import argparse,json
from pathlib import Path
import numpy as np
from hispid import Backend
from checkpoints import ROOT,select_record,restore
from physical import constraints,norms

p=argparse.ArgumentParser();p.add_argument('--library',required=True)
p.add_argument('--case',required=True);p.add_argument('--levels',required=True)
p.add_argument('--output',required=True);a=p.parse_args();b=Backend(a.library)
levels=[tuple(map(int,v.split(':'))) for v in a.levels.split(',')]
records=[];fixed=None
for n,np_ in levels:
    rec=select_record(a.case,n,np_);cfg,unknowns=restore(b,rec)
    centers=np.array([h.center for h in cfg.hole]);origin=centers.mean(axis=0)
    axis=centers[0]-centers[1];length=np.linalg.norm(axis);axis/=length
    helper=np.eye(3)[np.argmin(abs(axis))];y=np.cross(axis,helper);y/=np.linalg.norm(y);z=np.cross(axis,y)
    mass=sum(h.mass for h in cfg.hole)
    axial=[origin,origin+length*axis,origin-length*axis]
    for h in cfg.hole:
        if h.mass>0:
            for sign in (-1,1):axial.append(np.array(h.center)+sign*.6*mass*axis)
    x=np.array([point+radius*mass*(np.cos(phi)*y+np.sin(phi)*z)
                for point in axial for radius in (0,1e-6,1e-3,.01,.1) for phi in (.37,1.11)])
    if fixed is not None and not np.array_equal(fixed,x):raise ValueError('physical geometry differs between grids')
    fixed=x;sequence=[]
    with b.create_sampler(cfg) as s:
        s.set_unknowns(unknowns)
        for step in (.004,.002,.001):
            result=constraints(s.sample,x,step*mass)
            if np.any(result['attenuation']!=1):raise ValueError('axis control entered modified region')
            sequence.append(dict(step=step*mass,norms=norms(result),
                exact_axis=norms(result,np.linalg.norm((x-origin)-((x-origin)@axis)[:,None]*axis,axis=1)<1e-13*mass)))
    last=sequence[-1]['norms']
    accepted=rec['diagnostics']['status']==0 and all(last[q]<limit for q,limit in (('H_rms',1e-6),('M_rms',1e-6),('H_max',1e-4),('M_max',1e-4)))
    row=dict(resolution=rec['resolution'],step_checks=sequence,passed_strict=bool(accepted))
    records.append(row);print(json.dumps(row),flush=True)
report=dict(case=a.case,library_sha256=b.library_sha256(),parameterization=b.parameterization(),
    collocation_maps=b.parameterization_maps(),points=fixed.tolist(),records=records,
    passed=bool(len(records)>=3 and all(v['passed_strict'] for v in records[-3:])),
    continuous_exterior_supremum_verified=False)
(ROOT/a.output).write_text(json.dumps(report,indent=2)+'\n')
raise SystemExit(0 if report['passed'] else 1)
