"""Fully specified local inputs; distinctions from source cases are explicit."""
import numpy as np
from hispid import Hole

def moderate(backend,n=16,nphi=8):
    c=backend.config();c.n[:]=[n,n,nphi];c.tolerance=1e-14
    c.hole[0]=Hole(.6,(3,0,0),(.072,.054,.108),(.03,.06,.01))
    c.hole[1]=Hole(.4,(-3,0,0),(-.032,.048,.016),(-.02,-.07,.025))
    # Different spin axes, boost axes, unequal masses; |chi|=.3905/.3742.
    c.inner_min[:]=[.045,.03];c.inner_max[:]=[.09,.06]
    c.omega[:]=[.5,.5];c.attenuation_power=4;c.far_radius=40
    return c

def hs99uu(backend,n=32,nphi=16):
    c=backend.config();c.n[:]=[n,n,nphi];c.tolerance=1e-14
    c.hole[0]=Hole(.5,(6,0,0),(0,0,.2475));c.hole[1]=Hole(.5,(-6,0,0),(0,0,.2475))
    c.conformal_choice=0;c.omega[:]=[.2,.2];c.attenuation_power=4
    # Thesis chapter3 uses f only (no inner or far g); exact QI psi.
    c.inner_min[:]=[0,0];c.inner_max[:]=[0,0];c.inner_flatten=0;c.far_radius=0
    return c

def highboost(backend,n=32,nphi=8,momentum_over_rest_mass=2):
    c=backend.config();c.n[:]=[n,n,nphi];c.tolerance=1e-14
    v=momentum_over_rest_mass/np.sqrt(1+momentum_over_rest_mass**2)
    c.hole[0]=Hole(.5,(6,0,0),velocity=(-v,0,0));c.hole[1]=Hole(.5,(-6,0,0),velocity=(v,0,0))
    c.conformal_choice=0;c.omega[:]=[1,1];c.attenuation_power=4
    rh=.25/np.sqrt(1+momentum_over_rest_mass**2)
    c.inner_min[:]=[.2*rh,.2*rh];c.inner_max[:]=[.4*rh,.4*rh]
    c.inner_flatten=1;c.far_radius=40
    # Local v3 benchmark, NOT thesis table4.3: latter gives measured masses,
    # original discontinuous stuffing, and d/M_ADM=100--400.
    return c

def spin95(backend,n=32,nphi=16):
    c=hs99uu(backend,n,nphi)
    for hole in c.hole:hole.spin[:]=[0,0,.95*hole.mass**2]
    return c

def target_binary(backend,n=32,nphi=16,speed=.885,spin=.0,generic=False):
    """Revised local target: actual correction operators and a wider g window.

    The g-ball uses the contracted Kerr throat as a screen only. Actual
    enclosure is checked with the finder; it is not implied by these inputs.
    """
    if not 0<=speed<1 or not 0<=spin<1:raise ValueError('invalid target speed/spin')
    c=backend.config();c.n[:]=[n,n,nphi];c.tolerance=1e-14
    axes=[np.array([.2,.3,.4]),np.array([.4,-.2,.3])]
    boosts=[np.array([-1.,.23,.11]),np.array([1.,-.17,.19])]
    for h in range(2):
        spin_axis=axes[h]/np.linalg.norm(axes[h]) if generic else np.array([0,0,1.])
        boost_axis=boosts[h]/np.linalg.norm(boosts[h]) if generic else np.array([-1. if h==0 else 1.,0,0])
        c.hole[h]=Hole(.5,(6 if h==0 else -6,0,0),.25*spin*spin_axis,speed*boost_axis)
    c.conformal_choice=0;c.omega[:]=[1,1];c.attenuation_power=4
    radius=.25*np.sqrt(1-spin*spin)*np.sqrt(1-speed*speed)
    c.inner_min[:]=[.2*radius]*2;c.inner_max[:]=[.8*radius]*2
    c.inner_flatten=0;c.far_radius=0
    c.krylov_restart=80;c.max_krylov=2400;c.max_newton=24
    return c

def as_dict(config):
    return {name:[as_dict(h) for h in value] if name=='hole'
            else list(value) if hasattr(value,'_length_') else value
            for name,_ in config._fields_ for value in [getattr(config,name)]}
