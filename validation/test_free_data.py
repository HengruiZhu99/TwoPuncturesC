"""Focused oracle checks for physical metric/gdot free data (no solve sweep)."""
import argparse
import ctypes as C
import json
import sys
import tempfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Backend,Wave22,FreeData
from checkpoint_export import write_checkpoint,read_checkpoint
from kerr_wave_case import config,gauge

def oracle(w,x):
    # Independent spherical-basis implementation of the SpEC commit, not the
    # Cartesian TT projection or native jets. Points here are outside cutoff.
    r=np.linalg.norm(x,axis=1);rho=np.hypot(x[:,0],x[:,1]);phi=np.arctan2(x[:,1],x[:,0])
    cp=np.cos(phi);sp=np.sin(phi);mu=x[:,2]/r;st=rho/r
    et=np.array([mu*cp,mu*sp,-st]).T;ep=np.array([-sp,cp,np.zeros(len(x))]).T
    T=et[:,:,None]*et[:,None,:]-ep[:,:,None]*ep[:,None,:]
    U=et[:,:,None]*ep[:,None,:]+ep[:,:,None]*et[:,None,:]
    y=r-w.r0+w.propagation*w.time;L=2*np.pi*w.cycles/abs(w.omega)
    W=np.tanh(y)+np.tanh(L-y);Wp=1/np.cosh(y)**2-1/np.cosh(L-y)**2
    phase=2*phi-w.handedness*(w.omega*y+np.pi/4)
    u=np.sqrt(5/(64*np.pi))*(1+mu)**2;v=np.sqrt(5/(64*np.pi))*(1-mu)**2
    scale=np.sqrt(12)*w.amplitude/r
    a=scale*(u+v)*W*np.cos(phase);b=-scale*(u-v)*W*np.sin(phase)
    da=w.propagation*scale*(u+v)*(Wp*np.cos(phase)+w.handedness*w.omega*W*np.sin(phase))
    db=-w.propagation*scale*(u-v)*(Wp*np.sin(phase)-w.handedness*w.omega*W*np.cos(phase))
    return a[:,None,None]*T+b[:,None,None]*U,da[:,None,None]*T+db[:,None,None]*U

def run(b):
    rng=np.random.default_rng(4323);x=rng.normal(size=(80,3));x*=rng.uniform(4.1,25,(len(x),1))/np.linalg.norm(x,axis=1)[:,None]
    x=np.concatenate((x,[[0,0,10],[0,0,-10],[1e-14,0,10],[0,1e-14,-10]]))
    errors=[]
    for propagation in (-1,1):
        for handedness in (-1,1):
            w=Wave22(amplitude=.001,time=.37,propagation=propagation,handedness=handedness)
            a=b.wave22(w,x);ref=oracle(w,x)
            errors.append(max(np.max(abs(a[k]-ref[k])) for k in (0,1)))
    assert max(errors)<1e-17,errors
    w=Wave22(amplitude=.001);c=config(b)
    # Includes cutoff and axes, and a rotated solver frame (focus axis lab z).
    points=np.concatenate((x,[[0,0,.4],[2.5,.3,.2],[3.1,-.4,.5]]))
    dt=1e-5;wp=Wave22(**dict(w.as_dict(),time=dt));wm=Wave22(**dict(w.as_dict(),time=-dt))
    h,ht=b.wave22(w,points);fdtime=(b.wave22(wp,points)[0]-b.wave22(wm,points)[0])/(2*dt)
    time_error=float(np.max(abs(fdtime-ht)));assert time_error<2e-14,time_error
    with b.create_sampler(c) as old,b.create_sampler(c,wave=Wave22(amplitude=0)) as zero,b.create_sampler(c,wave=w) as s:
        a=old.sample(points);z=zero.sample(points);p=s.sample(points)
        assert all(np.array_equal(a[k],z[k]) for k in a)
        metric_error=float(np.max(abs(p['gamma'].reshape(-1,3,3)-a['gamma'].reshape(-1,3,3)-h)))
        assert metric_error<3e-14,metric_error
        # Reconstruct the ADM map independently using Cartesian differences.
        y=x[:8];hh,hh_t=b.wave22(w,y);alpha,beta=gauge(y);dh=[];db=[];eps=1e-4
        for d in range(3):
            shift=np.eye(3)[d]*eps
            dh.append((b.wave22(w,y+shift)[0]-b.wave22(w,y-shift)[0])/(2*eps))
            db.append((gauge(y+shift)[1]-gauge(y-shift)[1])/(2*eps))
        dh=np.array(dh);db=np.array(db)
        lie=np.einsum('nk,knij->nij',beta,dh)+np.einsum('nkj,ink->nij',hh,db)+np.einsum('nik,jnk->nij',hh,db)
        expected=-(hh_t-lie)/(2*alpha[:,None,None])
        actual=(s.sample(y)['Kij']-old.sample(y)['Kij']).reshape(-1,3,3)
        adm_error=float(np.max(abs(actual-expected)));assert adm_error<2e-12,adm_error
        # Native metric gradients use the same perturbed background as solve.
        grad=s.sample_with_derivatives(y)["dgamma"];eps=1e-4;gradient_error=0.
        for d in range(3):
            shift=np.eye(3)[d]*eps
            fd=(s.sample(y+shift)['gamma']-s.sample(y-shift)['gamma'])/(2*eps)
            gradient_error=max(gradient_error,float(np.max(abs(fd-grad[:,d].reshape(-1,9)))))
        assert gradient_error<1e-9,gradient_error
    tilted=config(b);tilted.hole[1].center[:]=(-1.,-.8,-1.2)
    with b.create_sampler(tilted,wave=w) as s,b.create_sampler(tilted) as old:
        a=old.sample(x[:8]);p=s.sample(x[:8]);h,_=b.wave22(w,x[:8])
        rotation_error=float(np.max(abs((p['gamma']-a['gamma']).reshape(-1,3,3)-h)))
        assert rotation_error<3e-14,rotation_error
        assert np.max(abs((p['Kij']-a['Kij']).reshape(-1,3,3)-expected))<2e-12
    # User g/gdot callback with analytic Cartesian jets, distinct from wave22.
    pol=np.polynomial.Polynomial.fromroots([2]*4+[5]*4)*1e-5
    def callback(x):
        r=np.linalg.norm(x);a=np.zeros((3,3,10));t=np.zeros_like(a)
        if 2<r<5:
            f=pol(r);fp=pol.deriv()(r);fpp=pol.deriv(2)(r);n=x/r
            Hess=(fpp-fp/r)*np.outer(n,n)+(fp/r)*np.eye(3)
            v=np.array([f,*(fp*n),Hess[0,0],Hess[0,1],Hess[0,2],Hess[1,1],Hess[1,2],Hess[2,2]])
            a[0,0]=v;a[1,1]=-v;t[0,1]=t[1,0]=.2*v
        return a,t
    provider=FreeData(callback)
    with b.create_sampler(c,free_data=provider) as s,b.create_sampler(c) as old:
        y=np.array([[3.,.2,.3]]);delta=callback(y[0])[0][:,:,0]
        assert np.max(abs((s.sample(y)['gamma']-old.sample(y)['gamma']).reshape(3,3)-delta))<2e-15
        try:write_checkpoint(Path(tempfile.gettempdir())/'should-not-write.hispid',s.config,s.unknowns(),b.library_sha256(),'diagnostic')
        except ValueError:pass
        else:raise AssertionError('callback data must not be silently dropped from export')
    # Reject data that would silently change the physical problem or geometry.
    bad=config(b);bad.hole[1].mass=.1
    for cfg,kwargs in ((bad,dict(wave=w)),(c,dict(wave=Wave22(cutoff_min=.1,cutoff_max=.3))),
                       (c,dict(wave=w,geometry='execution'))):
        try:b.create(cfg,**kwargs)
        except ValueError:pass
        else:raise AssertionError('unsupported free-data configuration accepted')
    for kind in ('nan','asymmetric','throat'):
        def invalid(x,kind=kind):
            a=np.zeros((3,3,10));t=np.zeros_like(a)
            if kind=='nan':a[0,0,0]=np.nan
            elif kind=='asymmetric':a[0,1,0]=1e-4
            else:a[0,0,0]=1e-4
            return a,t
        with b.create_sampler(c,free_data=FreeData(invalid)) as s:
            try:s.sample([[0,0,.35]])
            except ValueError:pass
            else:raise AssertionError('invalid callback accepted: '+kind)
    with b.create(c) as old,b.create(c,wave=Wave22(amplitude=0)) as zero:
        assert np.array_equal(old.residual(old.unknowns()),zero.residual(zero.unknowns()))
    # One tiny-grid algebra check, not a backend/performance sweep.
    with b.create(c,wave=w) as s:
        base=rng.normal(size=s.size)*1e-6;direction=rng.normal(size=s.size)*1e-6
        try:s.write_checkpoint(Path(tempfile.gettempdir())/'unsolved-wave.hispid')
        except ValueError:pass
        else:raise AssertionError('unsolved data must not be exported as a converged preliminary solve')
        j=s.jvp(base,direction);eps=1e-2
        fd=(s.residual(base+eps*direction)-s.residual(base-eps*direction))/(2*eps)
        jvp_error=float(np.max(abs(j-fd))/np.max(abs(j)));assert jvp_error<2e-8,jvp_error
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'wave.hispid';s.write_checkpoint(p,'diagnostic')
            cc,u,metadata=read_checkpoint(p);assert cc.wave22.as_dict()==w.as_dict()
            with b.create_sampler(cc) as loaded:
                loaded.set_unknowns(u)
                assert all(np.array_equal(s.sample(x[:2])[key],loaded.sample(x[:2])[key]) for key in ('gamma','Kij'))
    return dict(passed=True,oracle_absolute_linf=max(errors),fixed_coordinate_time_error=time_error,
                zero_amplitude_bitwise=True,metric_addition_error=metric_error,adm_map_error=adm_error,
                metric_gradient_error=gradient_error,arbitrary_frame_error=rotation_error,jvp_relative_error=jvp_error,checkpoint_roundtrip=True,
                library_sha256=b.library_sha256())

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path)
    a=p.parse_args();result=run(Backend(str(a.library.resolve())));text=json.dumps(result,indent=2);print(text)
    if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text+'\n')
