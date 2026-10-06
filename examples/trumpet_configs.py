"""Trumpet free data with explicitly enclosed candidate interior regularization.

Radii are seed-based screens, not measured enclosure certificates. The g=0
plateau plus inner_flatten=1 gives flat homogeneous correction equations at
both excluded foci; the physical trumpet seed itself is not flattened.
"""
import numpy as np
from configs import moderate

def trumpet_moderate(backend,n=24,nphi=8):
    c=moderate(backend,n,nphi);c.seed_family='trumpet_r0_m';c.n[:]=[n,2*n,nphi]
    c.conformal_choice=0;c.far_radius=0;c.inner_flatten=1
    c.tolerance=1e-12;c.max_newton=24;c.max_krylov=2400;c.krylov_restart=80
    for i,h in enumerate(c.hole):
        chi=np.linalg.norm(h.spin)/h.mass**2
        radius=h.mass*np.sqrt(1-chi**2)*np.sqrt(1-np.dot(h.velocity,h.velocity))
        c.inner_min[i]=.15*radius;c.inner_max[i]=.5*radius
    return c
