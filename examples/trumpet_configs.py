"""Trumpet binaries with explicit candidate interior regularization.

Radii screen for possible enclosure; only measured horizons can certify it.
The g=0 plateau plus inner_flatten=1 makes correction equations flat near
both foci. The physical trumpet seeds themselves are not flattened.
"""
import numpy as np
from hispid import Hole
from configs import moderate


def regularize_trumpet(config,n,nphi):
    c=config;c.seed_family='trumpet_r0_m';c.n[:]=[n,2*n,nphi]
    c.conformal_choice=0;c.far_radius=0;c.inner_flatten=1
    c.tolerance=1e-12;c.max_newton=24;c.max_krylov=2400;c.krylov_restart=80
    for i,h in enumerate(c.hole):
        if h.mass<=0:raise ValueError('binary requires two positive masses')
        chi=np.linalg.norm(h.spin)/h.mass**2;v2=np.dot(h.velocity,h.velocity)
        if not 0<=chi<1 or not 0<=v2<1:raise ValueError('subextremal spin and subluminal boost required')
        radius=h.mass*np.sqrt(1-chi**2)*np.sqrt(1-v2)
        c.inner_min[i]=.15*radius;c.inner_max[i]=.5*radius
    return c


def trumpet_moderate(backend,n=24,nphi=8):
    return regularize_trumpet(moderate(backend,n,nphi),n,nphi)


def trumpet_target(backend,case,n,nphi,separation=None,measured_mirr=None):
    """Separate spin/boost targets; boost separation needs measured calibration.

    measured_mirr is one component's measured irreducible mass. Its use sets
    the next trial's coordinate separation to50*Mirr; it does not certify the
    new trial's ratio, which must be measured again. Without a measurement,
    the boost trial starts at the exact isolated nonspinning value Mirr=.5.
    """
    if case not in ('spin99','gamma10'):raise ValueError('unknown separate target')
    if measured_mirr is not None:
        if case!='gamma10' or separation is not None:raise ValueError('mass calibration applies only to the boost separation')
        if not np.isfinite(measured_mirr) or measured_mirr<=0:raise ValueError('positive measured irreducible mass required')
        separation=50*measured_mirr
    if separation is None:separation=12. if case=='spin99' else 25.
    if not np.isfinite(separation) or separation<=0:raise ValueError('positive separation required')
    c=backend.config();c.omega[:]=[.5,.5];c.attenuation_power=4
    for i,sign in enumerate((1.,-1.)):
        c.hole[i]=Hole(.5,center=(sign*separation/2,0,0),
            spin=(0,0,.99*.5**2 if case=='spin99' else 0),
            velocity=(-sign*np.sqrt(.99) if case=='gamma10' else 0,0,0))
    return regularize_trumpet(c,n,nphi)
