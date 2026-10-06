"""Exact tensor-product interpolation of auxiliary grid unknowns.

Used only as a finer solve's initial guess; every level still solves its own
coupled equations and must pass the same residual and physical checks.
"""
import numpy as np

def chebyshev_evaluation_matrix(old,points):
    theta=np.pi*(np.arange(old)+.5)/old;x=-np.cos(theta)
    y=np.asarray(points,dtype=float)
    if y.ndim!=1 or not np.isfinite(y).all() or np.max(abs(y))>1+1e-14:
        raise ValueError('Chebyshev evaluation points must lie in[-1,1]')
    y=np.clip(y,-1,1)
    weights=(-1.)**np.arange(old)*np.sin(theta)
    out=np.empty((len(y),old))
    for i,point in enumerate(y):
        close=np.flatnonzero(np.abs(point-x)<2e-15)
        if len(close):out[i]=0;out[i,close[0]]=1
        else:
            row=weights/(point-x);out[i]=row/np.sum(row)
    return out

def chebyshev_matrix(old,new):
    return chebyshev_evaluation_matrix(old,-np.cos(np.pi*(np.arange(new)+.5)/new))

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

def prolong_modal(values,oldshape,newshape):
    """Regular modal P interpolation; phi slot is an orthonormal mode index."""
    a,b,p=oldshape;A,B,P=newshape
    if P<p:raise ValueError('modal prolongation must not discard Fourier modes')
    v=np.asarray(values).reshape(p,b,a,4)
    v=np.einsum('kjiv,ai->kjav',v,chebyshev_matrix(a,A),optimize=True)
    v=np.einsum('kjav,bj->kbav',v,chebyshev_matrix(b,B),optimize=True)
    return redistribute_modes(v,p,P).ravel()

def redistribute_modes(v,p,P):
    if p%2 or P%2 or P<p:raise ValueError('even Fourier grids without discarded modes required')
    out=np.zeros((P,*v.shape[1:]));half=p//2;newhalf=P//2
    for k in range(p):
        m=k if k<=half else k-half
        dest=m if k<=half else newhalf+m
        oldnormal=np.sqrt((1 if m in (0,half) else 2)/p)
        newnormal=np.sqrt((1 if m in (0,newhalf) else 2)/P)
        out[dest]=v[k]*oldnormal/newnormal
    return out

def validate_modal_shape(shape):
    if (len(shape)!=3 or any(isinstance(n,(bool,np.bool_)) or not isinstance(n,(int,np.integer)) for n in shape)
            or min(shape[:2])<2 or shape[2]<4 or shape[2]%2):
        raise ValueError('three integer grid dimensions with even Fourier dimension required')
    return tuple(shape)

def validate_modal_vector(values,shape):
    shape=validate_modal_shape(shape);values=np.asarray(values,dtype=float)
    if values.size!=4*int(np.prod(shape)) or not np.isfinite(values).all():
        raise ValueError('finite modal vector with the source grid coefficient count required')
    return values

def remap_modal(values,oldshape,newshape,oldmaps,newmaps):
    """Evaluate retained P directly at new physical t/eta coordinates.

    This is an initial-guess projection, never an accepted geometry migration.
    It avoids dividing physical fields by vanishing modal axis factors.
    """
    oldshape=validate_modal_shape(oldshape);newshape=validate_modal_shape(newshape)
    values=validate_modal_vector(values,oldshape)
    if len(oldmaps)!=2 or len(newmaps)!=2:raise ValueError('two map parameters required')
    for radial,angular in (oldmaps,newmaps):
        if not np.isfinite([radial,angular]).all() or not .001<=radial<=1 or not .1<=angular<=6:
            raise ValueError('unsupported collocation maps')
    if tuple(oldmaps)==tuple(newmaps):return prolong_modal(values,oldshape,newshape)
    a,b,p=oldshape;A,B,P=newshape
    if P<p:raise ValueError('remapping must not discard Fourier modes')
    radial_old,angular_old=oldmaps;radial_new,angular_new=newmaps
    sigma=.5*(1-np.cos(np.pi*(np.arange(A)+.5)/A))
    t=radial_new*sigma/(1-(1-radial_new)*sigma)
    radial_points=2*t/(radial_old+(1-radial_old)*t)-1
    zeta=-np.cos(np.pi*(np.arange(B)+.5)/B)
    eta=np.tanh(angular_new*zeta)/np.tanh(angular_new)
    angular_points=np.arctanh(np.tanh(angular_old)*eta)/angular_old
    v=np.asarray(values).reshape(p,b,a,4)
    v=np.einsum('kjiv,ai->kjav',v,chebyshev_evaluation_matrix(a,radial_points),optimize=True)
    v=np.einsum('kjav,bj->kbav',v,chebyshev_evaluation_matrix(b,angular_points),optimize=True)
    return redistribute_modes(v,p,P).ravel()

def for_backend(backend,values,oldshape,newshape):
    if backend.parameterization() in ('modal_P_C2prolate_v1','modal_P_C2prolate_mapped_v2') or backend.parameterization().startswith(('modal_P_C2prolate_map_v3_','modal_P_C4prolate_map_v4_')):return prolong_modal(values,oldshape,newshape)
    return prolong(values,oldshape,newshape)
