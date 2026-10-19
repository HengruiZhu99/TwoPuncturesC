"""Shared single-hole Liu configuration and independent stationary gauge."""
import numpy as np
from hispid import Hole


def config(b,n=8,chi=.9):
    c=b.config();c.hole[0]=Hole(1.,center=(0,0,0),spin=(0,0,chi))
    c.hole[1]=Hole(0.,center=(0,0,-2));c.seed_family='liu'
    c.n[:]=(n,2*n,16);c.far_radius=0;c.inner_min[:]=(.1,0);c.inner_max[:]=(.2,0)
    c.omega[:]=(0,0);c.memory_limit_mib=32768;c.tolerance=1e-12
    c.max_newton=12;c.max_krylov=1200;c.krylov_restart=100
    return c


def gauge(x,chi=.9):
    r=np.linalg.norm(x,axis=-1);c=(1+np.sqrt(1-chi*chi))/4;rb=r+2*c+c*c/r
    rm=chi*chi/(4*c);radial=(r-c)*np.sqrt((rb-rm)/r)
    sigma=rb*rb+chi*chi*(x[...,2]/r)**2
    AA=(rb*rb+chi*chi)**2-chi*chi*radial*radial*(1-(x[...,2]/r)**2)
    alpha=radial*np.sqrt(sigma/AA)
    cross=np.stack((-x[...,1],x[...,0],np.zeros_like(r)),axis=-1)
    return alpha,(-2*chi*rb/AA)[...,None]*cross

