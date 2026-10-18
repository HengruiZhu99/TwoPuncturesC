"""Independent analytic curvature and momentum controls for the FD6 verifier."""
import json,sys,argparse
from pathlib import Path
import numpy as np
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'validation'))
from physical import constraints,norms
p=argparse.ArgumentParser();p.add_argument('--points',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
# Rational Taylor moments through the formal accuracy orders.
x=(-3,-2,-1,1,2,3);d1=[F(q,60) for q in (-1,9,-45,45,-9,1)];d2=[F(q,180) for q in (2,-27,270,270,-27,2)]
for k in range(7):assert sum(w*F(z)**k for w,z in zip(d1,x))==int(k==1)
for k in range(8):assert sum(w*F(z)**k for w,z in zip(d2,x))-(F(49,18) if k==0 else 0)==2*int(k==2)
with np.load(a.points) as z:xyz=z['xyz'][:66];step=z['step'][:66]
def bl(x):
 psi=1+.3/np.linalg.norm(x-[3,0,0],axis=1)+.2/np.linalg.norm(x-[-3,0,0],axis=1)
 return dict(gamma=psi[:,None,None]**4*np.eye(3),Kij=np.zeros((len(x),3,3)),attenuation=np.ones(len(x)))
records={}
for label,order,factor in [('original_fd4',4,1),('declared_fd6',6,4)]:
 q=constraints(bl,xyz,factor*step,order=order);records[label]={name:norms(q,mask) for name,mask in [('near',np.arange(66)<48),('bulk',np.arange(66)>=48)]}
assert records['declared_fd6']['near']['H_rms']<1e-8
assert records['declared_fd6']['bulk']['H_rms']<1e-10
assert records['declared_fd6']['near']['M_rms']==0
# Flat gamma, Kij=f delta gives H=6f^2 and Mi=-2 partial_i f exactly.
xx=np.array([[.7,.4,.3],[.4,-.2,.1]])
def manufactured(x):
 f=.001*(x[:,0]**4+.5*x[:,1]**3-x[:,2]**5)
 return dict(gamma=np.broadcast_to(np.eye(3),(len(x),3,3)),Kij=f[:,None,None]*np.eye(3),attenuation=np.ones(len(x)))
q=constraints(manufactured,xx,.01,order=6);f=.001*(xx[:,0]**4+.5*xx[:,1]**3-xx[:,2]**5);M=-.002*np.column_stack([4*xx[:,0]**3,1.5*xx[:,1]**2,-5*xx[:,2]**4])
error=float(max(np.max(abs(q['H']-6*f*f)),np.max(abs(q['M']-M))));assert error<1e-13,error
a.output.write_text(json.dumps(dict(passed=True,calibration=records,manufactured_absolute_error=error),indent=2)+'\n');print(records['declared_fd6'])
