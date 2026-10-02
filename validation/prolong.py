"""Exact tensor-product interpolation of auxiliary grid unknowns.

Used only as a finer solve's initial guess; every level still solves its own
coupled equations and must pass the same residual and physical checks.
"""
import numpy as np

def chebyshev_matrix(old,new):
    theta=np.pi*(np.arange(old)+.5)/old;x=-np.cos(theta)
    y=-np.cos(np.pi*(np.arange(new)+.5)/new)
    weights=(-1.)**np.arange(old)*np.sin(theta)
    out=np.empty((new,old))
    for i,point in enumerate(y):
        close=np.flatnonzero(np.abs(point-x)<2e-15)
        if len(close):out[i]=0;out[i,close[0]]=1
        else:
            row=weights/(point-x);out[i]=row/np.sum(row)
    return out

def fourier_matrix(old,new):
    delta=2*np.pi*(np.arange(new)[:,None]/new-np.arange(old)[None,:]/old)
    out=np.ones((new,old))+np.cos(old/2*delta)
    for l in range(1,old//2):out+=2*np.cos(l*delta)
    return out/old

def prolong(values,oldshape,newshape):
    a,b,p=oldshape;A,B,P=newshape
    v=np.asarray(values).reshape(p,b,a,4)
    v=np.einsum('kjiv,ai->kjav',v,chebyshev_matrix(a,A),optimize=True)
    v=np.einsum('kjav,bj->kbav',v,chebyshev_matrix(b,B),optimize=True)
    v=np.einsum('kbav,ck->cbav',v,fourier_matrix(p,P),optimize=True)
    return np.ascontiguousarray(v).ravel()
