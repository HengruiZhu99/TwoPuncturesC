"""Manufactured curved scalar/vector operators with an independent FD observer.

The observer only reads conformal metric values. It independently constructs
connections and contravariant longitudinal tensors, then differences their
fluxes. Four distinct sinusoidal fields supply exact input Cartesian jets.
"""
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
from hispid import Backend,Hole
from physical import constraints

def run(library):
    b=Backend(str(library.resolve()));c=b.config();c.seed_family='trumpet_r0_m'
    c.hole[0]=Hole(1,spin=(.2,-.3,.7),velocity=(.25,.1,-.2));c.hole[1].mass=0
    c.omega[:]=[0,0];c.far_radius=0;c.inner_flatten=0;c.inner_min[:]=[0,0];c.inner_max[:]=[0,0]
    wave=np.array([[.7,-.2,.5],[.3,.6,-.2],[-.4,.5,.3],[.2,-.3,.8]])
    phase=np.arange(4)*.2
    def fields(x):
        t=x@wave.T+phase
        return np.sin(t),np.cos(t)[:,:,None]*wave[None],-np.sin(t)[:,:,None,None]*wave[None,:,:,None]*wave[None,:,None,:]
    def metric(x): return b.seed(c.hole[0],x,choice=c.conformal_choice,seed_family=c.seed_family)['conformal_metric'].reshape(-1,3,3)
    def derivative(fun,x,h):
        return np.stack([(fun(x-2*h*e)-8*fun(x-h*e)+8*fun(x+h*e)-fun(x+2*h*e))/(12*h) for e in np.eye(3)],axis=1)
    def geometry(x,h):
        g=metric(x);inv=np.linalg.inv(g);dg=derivative(metric,x,h)
        conn=np.empty((len(x),3,3,3))
        for k in range(3):
            for i in range(3):
                for j in range(3):conn[:,k,i,j]=.5*np.einsum('nl,nl->n',inv[:,k],dg[:,i,j]+dg[:,j,i]-dg[:,:,i,j])
        return inv,conn
    def longitudinal(x,h):
        inv,conn=geometry(x,h);v,dv,_=fields(x)
        D=np.transpose(dv[:,1:,:],(0,2,1))+np.einsum('njkl,nl->nkj',conn,v[:,1:])
        raised=inv@D
        return raised+raised.transpose(0,2,1)-(2/3)*inv*np.trace(D,axis1=1,axis2=2)[:,None,None]
    x=np.array([[1.2,.7,-.4],[2,-1,.8],[4,.3,2.]])
    v,dv,ddv=fields(x);jets=np.zeros((len(x),4,10));jets[:,:,0]=v;jets[:,:,1:4]=dv
    for k,(i,j) in enumerate([(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]):jets[:,:,4+k]=ddv[:,:,i,j]
    exact=np.array([b.operators(c,p,j) for p,j in zip(x,jets)])
    rows=[]
    for h in [.032,.016,.008]:
        inv,conn=geometry(x,h)
        lap=np.einsum('nij,nij->n',inv,ddv[:,0]-np.einsum('nkij,nk->nij',conn,dv[:,0]))
        L=longitudinal(x,h);dL=derivative(lambda y:longitudinal(y,h),x,h)
        div=np.einsum('njij->ni',dL)+np.einsum('nijk,nkj->ni',conn,L)+np.einsum('njjk,nik->ni',conn,L)
        def raw(y):return dict(gamma=metric(y),Kij=np.zeros((len(y),3,3)),attenuation=np.ones(len(y)))
        R=constraints(raw,x,h)['R']
        observed=np.column_stack([lap,div,R]);err=np.abs(observed-exact)/(1+np.abs(exact))
        rows.append(dict(step=h,component_scaled_linf=np.max(err,axis=0).tolist()))
    errors=np.array([row['component_scaled_linf'] for row in rows])
    return dict(kind='trumpet_manufactured_curved_operators',library_sha256=b.library_sha256(),
                rows=rows,passed=bool(np.all(errors[-1]<1e-7) and np.all(errors[1:]<errors[:-1])),
                relative_tolerance=1e-7,binary_acceptance=False,
                driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--library',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=run(a.library);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
