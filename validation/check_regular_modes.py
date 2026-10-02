"""Independent Cartesian scalar/vector oracle, including all axis segments.

The complex-distance expression uses physical Cartesian coordinates, without
the native compactification or derivative routines.
"""
import argparse,json
from pathlib import Path
import numpy as np
from hispid import Backend,Hole
from physical import constraints
from cartesian_modes import oracle

p=argparse.ArgumentParser();p.add_argument('--library',required=True);p.add_argument('--output',default='validation/regular_modes.json');a=p.parse_args()
b=Backend(a.library);c=b.config();n=64;np_=16;half=np_//2;c.n[:]=[n,n,np_];c.conformal_choice=0;c.inner_flatten=0
c.hole[0]=Hole(.6,(3,0,0));c.hole[1]=Hole(.4,(-3,0,0));c.omega[:]=[0,0];c.inner_max[:]=[0,0];c.far_radius=0
separation=3.;amp=.0001;zt=-np.cos(np.pi*(np.arange(n)+.5)/n);sigma=(1+zt)/2;maps=b.parameterization_maps()
lam=maps['radial_stretch'];kap=maps['angular_stretch'];t=lam*sigma/(1-(1-lam)*sigma);eta=np.tanh(kap*zt)/np.tanh(kap) if kap else zt
x=np.array([[xx,rr*np.cos(phi),rr*np.sin(phi)] for xx in (0.,1.2,5.4,-5.4) for rr in [0.,*[separation*10.**(-k) for k in range(2,11)],.1*separation,.6*separation,separation,2*separation] for phi in (.37,1.11)])
axis=np.linalg.norm(x[:,1:],axis=1)==0;ordinary=np.linalg.norm(x[:,1:],axis=1)>.5*separation;check=np.flatnonzero(axis|ordinary);records=[]
with b.create_sampler(c) as s:
 for m in (0,1,2,3,4,5,6,8):
  exponent=m if m<=4 else 3 if m%2 else 4
  for sine in (False,True):
   if sine and m in (0,half):continue
   mode=m if not sine else half+m;normal=np.sqrt((1 if m in (0,half) else 2)/np_)
   P=amp*((t[None,:]*(1-eta[:,None]**2))**((m-exponent)//2))/normal
   v,dv,ddv,seed,dseed=oracle(x,c.hole,separation,amp,m,sine)
   for component in (0,1):
    value=np.zeros((np_,n,n,4));value[mode,:,:,component]=P;s.set_unknowns(value.ravel());out=s.sample_with_derivatives(x)
    rec=dict(m=m,sine=sine,component='u' if component==0 else 'bx',maximum_correction_error=float(np.max(abs(out['correction'][:,component]-v))))
    if component==0:
     psi=seed+v;expected=4*psi[:,None,None,None]**3*(dseed+dv)[:,:,None,None]*np.eye(3)
     expectedH=-8*np.trace(ddv,axis1=1,axis2=2)/psi**5;expectedM=np.zeros((len(x),3));expectedK=np.zeros((len(x),3,3))
    else:
     psi=seed;expected=4*psi[:,None,None,None]**3*dseed[:,:,None,None]*np.eye(3)
     L=np.zeros((len(x),3,3));L[:,0,0]=4*dv[:,0]/3;L[:,1,1]=L[:,2,2]=-2*dv[:,0]/3
     L[:,0,1]=L[:,1,0]=dv[:,1];L[:,0,2]=L[:,2,0]=dv[:,2];expectedK=L/psi[:,None,None]**2
     expectedH=-((8/3)*dv[:,0]**2+2*dv[:,1]**2+2*dv[:,2]**2)/psi**12
     expectedM=ddv[:,0,:]/3;expectedM[:,0]+=np.trace(ddv,axis1=1,axis2=2);expectedM/=psi[:,None]**10
    rec['maximum_metric_gradient_error']=float(np.max(abs(out['dgamma']-expected)))
    rec['scaled_metric_gradient_error']=float(np.max(abs(out['dgamma']-expected)/(1+abs(expected))))
    rec['maximum_extrinsic_curvature_error']=float(np.max(abs(out['Kij'].reshape(-1,3,3)-expectedK)));rec['physical_step_checks']=[]
    for step in (.004,.002):
     actual=constraints(s.sample,x[check],step)
     rec['physical_step_checks'].append(dict(step=step,maximum_H_error=float(np.max(abs(actual['H']-expectedH[check]))),maximum_M_error=float(np.max(abs(actual['M']-expectedM[check])))))
    last=rec['physical_step_checks'][-1]
    rec['passed']=bool(rec['maximum_correction_error']<1e-12 and rec['scaled_metric_gradient_error']<1e-10 and rec['maximum_extrinsic_curvature_error']<1e-12 and last['maximum_H_error']<1e-8 and last['maximum_M_error']<1e-8)
    records.append(rec);print(rec,flush=True)
report=dict(library_sha256=b.library_sha256(),parameterization=b.parameterization(),collocation_maps=maps,resolution=list(c.n),points=x.tolist(),physical_check_indices=check.tolist(),amplitude=amp,records=records,passed=all(v['passed'] for v in records))
Path(a.output).write_text(json.dumps(report,indent=2)+'\n');raise SystemExit(0 if report['passed'] else 1)
