"""Cartesian distance-mode value/gradient/Hessian oracle, independent of native maps."""
import numpy as np

def oracle(x,holes,separation,amp,m,sine):
    delta=[x-np.array(h.center) for h in holes]
    radii=[np.linalg.norm(d,axis=1) for d in delta]
    unit=[d/r[:,None] for d,r in zip(delta,radii)]
    D=(sum(radii)+2*separation)/2;dD=.5*sum(unit)
    ddD=.5*sum((np.eye(3)-u[:,:,None]*u[:,None,:])/r[:,None,None] for u,r in zip(unit,radii))
    w=x[:,1]+1j*x[:,2];wi=np.array([0,1,1j]);z=w**m
    pref=-4*separation*amp;take=np.imag if sine else np.real
    value=pref*take(z/D**(m+1))
    grad=-(m+1)*z[:,None]*dD/D[:,None]**(m+2)
    hess=((m+1)*(m+2)*z[:,None,None]*dD[:,:,None]*dD[:,None,:]/D[:,None,None]**(m+3)
          -(m+1)*z[:,None,None]*ddD/D[:,None,None]**(m+2))
    if m:
        grad+=m*w[:,None]**(m-1)*wi/D[:,None]**(m+1)
        hess-=m*(m+1)*w[:,None,None]**(m-1)*(wi[None,:,None]*dD[:,None,:]+dD[:,:,None]*wi[None,None,:])/D[:,None,None]**(m+2)
    if m>=2:hess+=m*(m-1)*w[:,None,None]**(m-2)*wi[None,:,None]*wi[None,None,:]/D[:,None,None]**(m+1)
    seed=1+sum(h.mass/(2*r) for h,r in zip(holes,radii))
    dseed=-sum(h.mass*d/(2*r[:,None]**3) for h,d,r in zip(holes,delta,radii))
    return value,pref*take(grad),pref*take(hess),seed,dseed

