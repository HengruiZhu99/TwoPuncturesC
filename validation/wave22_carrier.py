"""Independent analytic complex radial carrier and its first two derivatives.

The native physical tensor is tested against these formulas by the profile
oracle. Spectral resolution is assessed in physical radial rays by
wave22_physical_radial_audit.py, including both prolate coordinates.
"""
import numpy as np


def carrier(w,r):
    y=r-w.r0+w.propagation*w.time;end=2*np.pi*w.cycles/abs(w.omega)-y
    # Stable tanh window and derivatives; no overflow at compactified infinity.
    def sech2(x):
        e=np.exp(-2*abs(x));return 4*e/(1+e)**2
    W=np.tanh(y)+np.tanh(end)
    wp=sech2(y)-sech2(end);wpp=-2*np.tanh(y)*sech2(y)-2*np.tanh(end)*sech2(end)
    k=w.handedness*w.omega;phase=np.exp(-1j*(k*y+w.handedness*np.pi/4))
    return np.stack((W*phase,(wp-1j*k*W)*phase,(wpp-2j*k*wp-k*k*W)*phase),axis=-1)
