"""Independent Cartesian checks of raw physical gamma_ij and K_ij.

No native conformal derivatives, constraints, connections, or curvature reused.
Fourth-order Cartesian differences, including tensor-product mixed derivatives.
"""
from __future__ import annotations
import numpy as np

def constraints(sample,xyz,step):
    x=np.array(xyz,dtype=float).reshape(-1,3);n=len(x)
    steps=np.broadcast_to(np.asarray(step,dtype=float),(n,))
    if np.any(steps<=0) or not np.all(np.isfinite(steps)):raise ValueError('finite positive verifier steps required')
    offsets=[np.zeros(3)]
    axial={};mixed={}
    for d in range(3):
        axial[d]=[]
        for a in (-2,-1,1,2):
            dx=np.zeros(3);dx[d]=a;axial[d].append(len(offsets));offsets.append(dx)
    for i in range(3):
        for j in range(i+1,3):
            mixed[i,j]=[]
            for a in (-2,-1,1,2):
                for b in (-2,-1,1,2):
                    dx=np.zeros(3);dx[i]=a;dx[j]=b
                    mixed[i,j].append(len(offsets));offsets.append(dx)
    pts=x[None,:,:]+steps[None,:,None]*np.array(offsets)[:,None,:]
    values=sample(pts.reshape(-1,3))
    g=values['gamma'].reshape(-1,n,3,3);K=values['Kij'].reshape(-1,n,3,3)
    d1=np.array([1,-8,8,-1])/12
    d2=np.array([-1,16,16,-1])/12
    dg=np.zeros((n,3,3,3));ddg=np.zeros((n,3,3,3,3));dk=np.zeros_like(dg)
    for d in range(3):
        dg[:,d]=np.einsum('a,anij->nij',d1,g[axial[d]])/steps[:,None,None]
        dk[:,d]=np.einsum('a,anij->nij',d1,K[axial[d]])/steps[:,None,None]
        ddg[:,d,d]=(np.einsum('a,anij->nij',d2,g[axial[d]])-2.5*g[0])/(steps[:,None,None]**2)
    for (i,j),ix in mixed.items():
        v=np.einsum('a,anij->nij',np.outer(d1,d1).reshape(-1),g[ix])/(steps[:,None,None]**2)
        ddg[:,i,j]=ddg[:,j,i]=v
    inv=np.linalg.inv(g[0]);dinv=-np.einsum('nik,ndkl,nlj->ndij',inv,dg,inv)
    conn=np.zeros((n,3,3,3));dc=np.zeros((n,3,3,3,3))
    for k in range(3):
        for i in range(3):
            for j in range(3):
                q=dg[:,i,j,:]+dg[:,j,i,:]-dg[:,:,i,j]
                conn[:,k,i,j]=.5*np.einsum('nl,nl->n',inv[:,k,:],q)
                for d in range(3):
                    dq=ddg[:,d,i,j,:]+ddg[:,d,j,i,:]-ddg[:,d,:,i,j]
                    dc[:,d,k,i,j]=.5*(np.einsum('nl,nl->n',dinv[:,d,k,:],q)+np.einsum('nl,nl->n',inv[:,k,:],dq))
    ricci=np.zeros((n,3,3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ricci[:,i,j]+=dc[:,k,k,i,j]-dc[:,j,k,i,k]
                for l in range(3):
                    ricci[:,i,j]+=conn[:,k,i,j]*conn[:,l,k,l]-conn[:,l,i,k]*conn[:,k,j,l]
    R=np.einsum('nij,nij->n',inv,ricci)
    tr=np.einsum('nij,nij->n',inv,K[0]);kmix=inv@K[0]
    K2=np.einsum('nij,nji->n',kmix,kmix)
    identity=np.eye(3);B=kmix-tr[:,None,None]*identity
    dtrace=np.einsum('ndij,nij->nd',dinv,K[0])+np.einsum('nij,ndij->nd',inv,dk)
    dB=np.einsum('ndij,njk->ndik',dinv,K[0])+np.einsum('nij,ndjk->ndik',inv,dk)-dtrace[:,:,None,None]*identity
    Mlow=np.zeros((n,3))
    for i in range(3):
        for j in range(3):
            Mlow[:,i]+=dB[:,j,j,i]
            for k in range(3):
                Mlow[:,i]+=conn[:,j,j,k]*B[:,k,i]-conn[:,k,j,i]*B[:,j,k]
    M=np.einsum('nij,nj->ni',inv,Mlow)
    Mnorm=np.sqrt(np.maximum(0,np.einsum('ni,nij,nj->n',Mlow,inv,Mlow)))
    H=R+tr*tr-K2
    divKlow=Mlow+dtrace
    divKnorm=np.sqrt(np.maximum(0,np.einsum('ni,nij,nj->n',divKlow,inv,divKlow)))
    gradKnorm=np.sqrt(np.maximum(0,np.einsum('ni,nij,nj->n',dtrace,inv,dtrace)))
    return dict(H=H,M=M,Mnorm=Mnorm,R=R,trace=tr,K2=K2,
                normalized_H=np.abs(H)/np.maximum(np.abs(R)+tr*tr+K2,1e-8),
                normalized_M=Mnorm/np.maximum(divKnorm+gradKnorm,1e-8),
                attenuation=values['attenuation'].reshape(-1,n)[0],
                stencil_attenuation_min=np.min(values['attenuation'].reshape(-1,n),axis=0),
                stencil_attenuation_all_one=np.all(values['attenuation'].reshape(-1,n)==1,axis=0),
                min_metric_eigenvalue=np.linalg.eigvalsh(g[0])[:,0])

def norms(r,mask=None):
    if mask is None:mask=np.ones(len(r['H']),dtype=bool)
    n=int(np.sum(mask))
    if not n:return {'count':0}
    return dict(count=n,H_rms=float(np.sqrt(np.mean(r['H'][mask]**2))),
                H_max=float(np.max(np.abs(r['H'][mask]))),
                M_rms=float(np.sqrt(np.mean(r['Mnorm'][mask]**2))),
                M_max=float(np.max(r['Mnorm'][mask])),
                M_components_rms=np.sqrt(np.mean(r['M'][mask]**2,axis=0)).tolist(),
                M_components_max=np.max(np.abs(r['M'][mask]),axis=0).tolist(),
                normalized_H_max=float(np.max(r['normalized_H'][mask])),
                normalized_M_max=float(np.max(r['normalized_M'][mask])))

def charges(sample,radius,center=(0,0,0),ntheta=20,nphi=40,polar_frame=None):
    """Independent finite-radius ADM integrals of physical fields, flat measure."""
    mu,w=np.polynomial.legendre.leggauss(ntheta);phi=2*np.pi*(np.arange(nphi)+.5)/nphi
    mu,phi=np.meshgrid(mu,phi,indexing='ij');normal=np.stack([np.sqrt(1-mu*mu)*np.cos(phi),np.sqrt(1-mu*mu)*np.sin(phi),mu],axis=-1).reshape(-1,3)
    if polar_frame is not None:
        frame=np.asarray(polar_frame,dtype=float)
        if frame.shape!=(3,3) or not np.allclose(frame.T@frame,np.eye(3),rtol=0,atol=1e-12) or abs(np.linalg.det(frame)-1)>1e-12:
            raise ValueError('polar_frame must be a proper orthonormal frame')
        # The first frame column is the polar axis; the other two span phi.
        normal=normal[:,[2,0,1]]@frame.T
    weight=np.repeat(w,nphi)*2*np.pi/nphi*radius**2
    x=np.asarray(center)+radius*normal;v=sample(x);g=v['gamma'].reshape(-1,3,3);K=v['Kij'].reshape(-1,3,3)
    inv=np.linalg.inv(g);tr=np.einsum('nij,nij->n',inv,K)
    h=radius*1e-4;dg=np.empty((len(x),3,3,3))
    for d in range(3):
        shifts=np.zeros((4,3));shifts[:,d]=np.array([-2,-1,1,2])*h
        gs=sample((x[None,:,:]+shifts[:,None,:]).reshape(-1,3))['gamma'].reshape(4,-1,3,3)
        dg[:,d]=np.einsum('a,anij->nij',np.array([1,-8,8,-1])/12,gs)/h
    flux=np.zeros((len(x),3))
    for i in range(3):
        for j in range(3):flux[:,i]+=dg[:,j,i,j]-dg[:,i,j,j]
    E=np.sum(weight*np.einsum('ni,ni->n',normal,flux))/(16*np.pi)
    pflux=np.einsum('nij,nj->ni',K-tr[:,None,None]*g,normal)
    P=np.einsum('n,ni->i',weight,pflux)/(8*np.pi)
    J=np.einsum('n,ni->i',weight,np.cross(radius*normal,pflux))/(8*np.pi)
    return np.r_[E,P,J]

def extrapolate(radii,charges_by_radius,degree=2):
    x=1/np.asarray(radii,dtype=float)
    return np.polynomial.polynomial.polyfit(x,np.asarray(charges_by_radius),degree)[0]
