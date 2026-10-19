"""Audit spectral carrier resolution along physical single-hole radial rays.

Physical radius changes both prolate coordinates. This two-dimensional test
therefore complements the fixed-polar radial-line diagnostic. It interpolates
the complex carrier at actual mapped Chebyshev nodes and differentiates along
fixed physical polar angle, against an independent analytic radial formula.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.fft import dct
from numpy.polynomial.chebyshev import chebval2d,chebder
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Wave22
from wave22_carrier import carrier


def audit(w,nr,ntheta,separation,stretch,points=513):
    b=separation/2;kappa=2.;T=np.tanh(kappa);hi=w.r0+2*np.pi*w.cycles/abs(w.omega)+4;lo=w.r0-4
    radii=np.unique(np.r_[np.linspace(lo,hi,points),[separation+e for e in (-.01,-1e-4,1e-4,.01) if lo<separation+e<hi]])
    rr,mu=np.meshgrid(radii,np.linspace(-1,1,9),indexing='ij');rr=rr.ravel();mu=mu.ravel()
    other=np.hypot(rr+separation*mu,separation*np.sqrt(1-mu*mu));mask=other>1e-8
    rr=rr[mask];mu=mu[mask];other=other[mask]
    C=(rr+other)/separation;eta=(other-rr)/separation;tt=(C-1)/(C+1);den=stretch+(1-stretch)*tt
    x=2*tt/den-1;zeta=np.arctanh(np.clip(eta*T,-T,T))/kappa
    dother=(rr+separation*mu)/other;Cr=(1+dother)/separation;er=(dother-1)/separation
    xr=2*stretch/den**2*2*Cr/(C+1)**2;zr=T*er/(kappa*(1-eta*eta*T*T))
    sigma=(1-np.cos(np.pi*(np.arange(nr)+.5)/nr))/2;t=stretch*sigma/(1-(1-stretch)*sigma)
    z=-np.cos(np.pi*(np.arange(ntheta)+.5)/ntheta);e=np.tanh(kappa*z)/T
    rnodes=b*((1+t[:,None])/(1-t[:,None])-e[None,:]);q=carrier(w,rnodes)[:,:,0]
    c=dct(dct(q[::-1,::-1],axis=0)/nr,axis=1)/ntheta;c[0,:]*=.5;c[:,0]*=.5
    observed=chebval2d(x,zeta,c)
    dr=chebval2d(x,zeta,chebder(c,axis=0))*xr+chebval2d(x,zeta,chebder(c,axis=1))*zr
    other_rr=separation**2*(1-mu*mu)/other**3
    Crr=other_rr/separation;err=other_rr/separation
    tr=2*Cr/(C+1)**2;trr=2*Crr/(C+1)**2-4*Cr*Cr/(C+1)**3
    xrr=2*stretch*(trr/den**2-2*(1-stretch)*tr*tr/den**3)
    zrr=T/kappa*(err/(1-eta*eta*T*T)+2*eta*T*T*er*er/(1-eta*eta*T*T)**2)
    qx=chebder(c,axis=0);qz=chebder(c,axis=1)
    dr2=chebval2d(x,zeta,chebder(qx,axis=0))*xr*xr+2*chebval2d(x,zeta,chebder(qx,axis=1))*xr*zr+chebval2d(x,zeta,chebder(qz,axis=1))*zr*zr+chebval2d(x,zeta,qx)*xrr+chebval2d(x,zeta,qz)*zrr
    ref=carrier(w,rr);errors=np.array([np.max(abs(observed-ref[:,0]))/2,np.max(abs(dr-ref[:,1]))/np.max(abs(ref[:,1])),np.max(abs(dr2-ref[:,2]))/np.max(abs(ref[:,2]))])
    active=abs(ref[:,0])>.2
    phase=float(np.max(abs(np.angle(observed[active]*np.conj(ref[active,0])))))
    return dict(n=[nr,ntheta,16],focus_separation=separation,radial_stretch=stretch,query_count=len(rr),excluded_degenerate_points=int(np.sum(~mask)),
                relative_linf=dict(value=float(errors[0]),dr=float(errors[1])),
                unconstrained_interpolant_dr2_relative_linf=float(errors[2]),phase_max_radians_above_ten_percent=phase,
                passed=bool(np.max(errors[:2])<1e-3 and phase<1e-3))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--grids',nargs='+',default=['160,320,32,.2','768,64,32,.2','384,128,16,1','480,128,16,1','576,128,16,1'])
    a=p.parse_args();d=json.loads(a.study.read_text());w=Wave22(**d['wave22']);rows=[]
    for grid in a.grids:
        nr,nb,sep,l=map(float,grid.split(','));row=audit(w,int(nr),int(nb),sep,l);rows.append(row);print(row,flush=True)
    result=dict(wave22=w.as_dict(),criteria=dict(value_relative_linf=1e-3,dr_relative_linf=1e-3,core_phase_radians=1e-3),
        definition='Complex carrier interpolated on the actual two-dimensional prolate grid; physical radial derivatives include both map coordinates; nine physical polar rays,513 radii plus one-sided probes around the inactive focus',
        limitations='Carrier value/first-derivative test, not full tensor-coefficient or constraint error bound; excludes only the coordinate-degenerate focus itself and probes both sides. The unconstrained scalar interpolant second derivative is diagnostic only: its tensor-product polynomial does not impose Cartesian focus regularity, while native free-data derivatives are analytic. It is not a bound on native tensor second derivatives; independent physical constraints are required',
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),study_sha256=hashlib.sha256(a.study.read_bytes()).hexdigest(),rows=rows,passed=bool(rows and rows[-1]['passed']))
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    if not result['passed']:raise SystemExit(1)
