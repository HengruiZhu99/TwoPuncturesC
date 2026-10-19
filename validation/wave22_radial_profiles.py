"""Dense physical wave profiles from immutable checkpoints; no new solves.

The seed is analytically sampled, so qualification also needs an independent
collocation carrier-resolution check. All differences here subtract Kerr and
are normalized to the wave, never to the much larger background metric.
"""
import argparse,copy,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Backend,Wave22
from kerr_wave_case import config
from checkpoint_export import read_checkpoint
from kerr_wave_case import gauge
from parallel_observer import ParallelObserver
from physical import constraints,norms
from wave22_carrier import carrier


def profile(s,radii,mu=(-.6,0.,.6),nphi=16,threads=16):
    mu=np.asarray(mu);phi=2*np.pi*np.arange(nphi)/nphi
    mm,pp=np.meshgrid(mu,phi,indexing='ij');sn=np.sqrt(1-mm*mm)
    normal=np.stack((sn*np.cos(pp),sn*np.sin(pp),mm),axis=-1).reshape(-1,3)
    et=np.stack((mm*np.cos(pp),mm*np.sin(pp),-sn),axis=-1).reshape(-1,3)
    ep=np.stack((-np.sin(pp),np.cos(pp),np.zeros_like(pp)),axis=-1).reshape(-1,3)
    x=(radii[:,None,None]*normal[None,:,:]).reshape(-1,3)
    n=np.tile(normal,(len(radii),1));et=np.tile(et,(len(radii),1));ep=np.tile(ep,(len(radii),1))
    chi=s.config.hole[0].spin[2]
    with ParallelObserver(s,threads) as sample:q=sample.with_derivatives(x)
    g=q['gamma'].reshape(-1,3,3);K=q['Kij'].reshape(-1,3,3);dg=q['dgamma']
    alpha,beta=gauge(x,chi);eps=1e-4
    db=np.array([(gauge(x+eps*np.eye(3)[d],chi)[1]-gauge(x-eps*np.eye(3)[d],chi)[1])/(2*eps) for d in range(3)])
    dt=np.einsum('nk,nkij->nij',beta,dg)+np.einsum('nkj,ink->nij',g,db)+np.einsum('nik,jnk->nij',g,db)-2*alpha[:,None,None]*K
    dr=np.einsum('nk,nkij->nij',n,dg)
    out={}
    for name,T in (('metric',g),('dt_metric',dt),('dr_metric',dr)):
        pol=.5*(np.einsum('ni,nij,nj->n',et,T,et)-np.einsum('ni,nij,nj->n',ep,T,ep))-1j*np.einsum('ni,nij,nj->n',et,T,ep)
        angular=pol.reshape(len(radii),len(mu),nphi)
        out[name]=np.stack([np.mean(angular*np.exp(-1j*m*phi)[None,None,:],axis=-1) for m in (2,-2)],axis=-1)
    return out


def derived(d,r):
    scale=r.reshape((len(r),)+(1,)*(d['metric'].ndim-1))
    h=scale*d['metric'];dt=scale*d['dt_metric'];dr=d['metric']+scale*d['dr_metric']
    return dict(h=h,dt_h=dt,dr_h=dr,ingoing=(dt+dr)/2,outgoing=(dt-dr)/2)


def compare(a,b,r):
    out={}
    for name in ('h','dt_h','dr_h','ingoing'):
        peak=np.max(abs(b[name]),axis=0);err=a[name]-b[name]
        out[name]=dict(relative_linf=float(np.max(np.max(abs(err),axis=0)/peak)),relative_l2=float(np.max(np.linalg.norm(err,axis=0)/np.linalg.norm(b[name],axis=0))))
    active=abs(b['h'])>.1*np.max(abs(b['h']),axis=0)
    out['phase_max_radians']=float(np.max(abs(np.angle(a['h'][active]*np.conj(b['h'][active])))))
    out['amplitude_peak_normalized_linf']=float(np.max(np.max(abs(abs(a['h'])-abs(b['h'])),axis=0)/np.max(abs(b['h']),axis=0)))
    out['outgoing_difference_over_ingoing_peak']=float(np.max(np.max(abs(a['outgoing']-b['outgoing']),axis=0)/np.max(abs(b['ingoing']),axis=0)))
    out['passed']=bool(all(out[k]['relative_linf']<1e-3 for k in ('h','dt_h','dr_h','ingoing')) and out['phase_max_radians']<1e-3 and out['outgoing_difference_over_ingoing_peak']<1e-3)
    return out


def self_test(b):
    c=config(b,chi=.9999);w=Wave22(amplitude=.001,omega=.9760972454792889,cycles=4.,r0=12.,cutoff_min=2.,cutoff_max=4.)
    radii=np.linspace(9,41,9);mu=np.array([-.6,0,.6]);Y=np.sqrt(5/(64*np.pi))*np.stack(((1+mu)**2,(1-mu)**2),axis=-1)
    with b.create_sampler(c) as seed,b.create_sampler(c,wave=w) as trial:
        baseline=profile(seed,radii,threads=4);p=profile(trial,radii,threads=4)
        actual={k:p[k]-baseline[k] for k in p}
    q=carrier(w,radii);q=np.stack((q,np.conj(q)),axis=-1);factor=np.sqrt(12)*w.amplitude*Y[None,:,:]/radii[:,None,None]
    expected=dict(metric=factor*q[:,0,None,:],dt_metric=factor*q[:,1,None,:],dr_metric=factor*(q[:,1,None,:]-q[:,0,None,:]/radii[:,None,None]))
    errors={key:float(np.max(abs(actual[key]-expected[key]))) for key in actual}
    assert max(errors.values())<1e-14,errors
    c.wave22=w;c0=copy.copy(c);del c0.wave22
    assert hasattr(c,'wave22') and not hasattr(c0,'wave22')
    return dict(passed=True,analytic_mode_linf=errors,library_sha256=b.library_sha256(),driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',type=Path,required=True)
    p.add_argument('--checkpoints',type=Path,nargs='+');p.add_argument('--output',type=Path,required=True)
    p.add_argument('--packet-constraints',action='store_true');p.add_argument('--self-test',action='store_true');p.add_argument('--points',type=int,default=257);p.add_argument('--threads',type=int,default=16);a=p.parse_args()
    b=Backend(str(a.library.resolve()))
    if a.self_test:
        result=self_test(b);a.output.write_text(json.dumps(result,indent=2)+'\n');print(result);return
    if not a.checkpoints:raise ValueError('checkpoints required')
    if a.points<65:raise ValueError('at least 65 radial observers required')
    a.output.mkdir(parents=True,exist_ok=False);rows=[];profiles=[];archive={}
    source_wave=None;packet_summary=None;source_case=None
    for path in a.checkpoints:
        c,u,meta=read_checkpoint(path)
        if meta['parameterization']!=b.parameterization():raise ValueError('checkpoint/basis mismatch')
        if meta['acceptance']=='diagnostic':raise ValueError('nonconverged diagnostic checkpoint cannot qualify the solved wave')
        if source_wave is not None and c.wave22.as_dict()!=source_wave:raise ValueError('different wave data')
        source_wave=c.wave22.as_dict();w=c.wave22
        case=dict(conformal_choice=c.conformal_choice,inner_flatten=c.inner_flatten,attenuation_power=c.attenuation_power,holes=[[h.mass,*h.center,*h.spin,*h.velocity] for h in c.hole],inner_min=list(c.inner_min),inner_max=list(c.inner_max),omega=list(c.omega),far_radius=c.far_radius,tolerance=c.tolerance,seed_family=c.seed_family)
        if source_case is not None and case!=source_case:raise ValueError('physical case or stopping tolerance changed')
        source_case=case
        if c.seed_family!='liu' or c.hole[0].mass!=1 or c.hole[1].mass!=0 or any(c.hole[0].velocity) or any(c.hole[0].center) or any(c.hole[0].spin[:2]) or not 0<c.hole[0].spin[2]<1:raise ValueError('centered mass-one stationary Liu hole with positive z spin required')
        if w.amplitude==0:raise ValueError('nonzero packet required for wave-normalized comparisons')
        radii=np.linspace(w.r0-4,w.r0+2*np.pi*w.cycles/abs(w.omega)+4,a.points)
        c0=copy.copy(c);del c0.wave22
        c0.n[:]=(8,8,8)  # exact zero-correction seed is independent of spectral size
        with b.create_sampler(c0) as seed,b.create_sampler(c) as solved:
            solved.set_unknowns(u);baseline=profile(seed,radii,threads=a.threads);raw=profile(solved,radii,threads=a.threads)
            if a.packet_constraints and path==a.checkpoints[-1]:
                packet_r=np.linspace(w.r0-2,w.r0+2*np.pi*w.cycles/abs(w.omega)+2,97)
                packet_mu=np.array([-1.,-.95,0.,.95,1.]);sn=np.sqrt(1-packet_mu*packet_mu)
                directions=np.stack((sn*np.cos(np.pi/8),sn*np.sin(np.pi/8),packet_mu),axis=-1)
                xyz=(packet_r[:,None,None]*directions[None,:,:]).reshape(-1,3);step=np.full(len(xyz),.005)
                with ParallelObserver(solved,a.threads) as sm:after=constraints(sm,xyz,step,6)
                with ParallelObserver(seed,a.threads) as sm:floor=norms(constraints(sm,xyz,step,6))
                trial_config=copy.copy(c);trial_config.n[:]=(8,8,8)
                with b.create_sampler(trial_config) as trial,ParallelObserver(trial,a.threads) as sm:before=norms(constraints(sm,xyz,step,6))
                actual=norms(after);packet_summary=dict(points=len(xyz),mu=packet_mu.tolist(),order=6,steps=.005,exact_seed_floor=floor,before=before,after=actual,
                    physical_maximum_bound=1e-8,absolute_bound_passed=bool(actual['H_max']<1e-8 and actual['M_max']<1e-8),
                    after_by_mu=[norms(after,np.arange(len(xyz))%len(packet_mu)==j) for j in range(len(packet_mu))],
                    projection_reduces_both=bool(actual['H_rms']<before['H_rms'] and actual['M_rms']<before['M_rms']))
                archive.update(packet_xyz=xyz,packet_H=after['H'],packet_Mnorm=after['Mnorm'])
        delta={key:raw[key]-baseline[key] for key in raw};v=derived(delta,radii);profiles.append(v)
        key='x'.join(map(str,c.n));archive.update({key+'_'+name:value for name,value in v.items()})
        row=dict(n=list(c.n),checkpoint=meta,ingoing_to_outgoing_note='coordinate metric characteristics, no curvature energy flux',
                 outgoing_over_ingoing_l2_by_mu=(np.linalg.norm(v['outgoing'],axis=0)/np.linalg.norm(v['ingoing'],axis=0)).tolist())
        if len(profiles)>1:row['previous_difference']=compare(profiles[-2],profiles[-1],radii)
        rows.append(row);print(key,row.get('previous_difference',{}),flush=True)
    decreasing=bool(len(rows)>=3 and all(rows[-1]['previous_difference'][key]['relative_linf']<rows[-2]['previous_difference'][key]['relative_linf'] for key in ('h','dt_h','dr_h','ingoing')))
    # Compare every level to the final profile without discarding earlier failures.
    for row,v in zip(rows,profiles):row['difference_to_finest']=compare(v,profiles[-1],radii)
    np.savez_compressed(a.output/'profiles.npz',radius=radii,mu=np.array([-.6,0,.6]),azimuthal_modes=np.array([2,-2]),**archive)
    result=dict(library_sha256=b.library_sha256(),driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        wave22=source_wave,points=a.points,radial_range=[float(radii[0]),float(radii[-1])],mu=[-.6,0,.6],azimuthal_modes=[2,-2],nphi=16,threads=a.threads,
        criteria=dict(wave_peak_normalized_linf=1e-3,phase_radians=1e-3),rows=rows,
        note='Kerr-subtracted m=+/-2 Fourier rings; relative errors normalized per ring to the finest wave peak, phase tested above10% peak; this profile check alone cannot exclude collocation aliasing',
        packet_constraints=packet_summary,profile_differences_decrease=decreasing,passed=bool(decreasing and rows[-1]['previous_difference']['passed'] and (packet_summary is None or (packet_summary['projection_reduces_both'] and packet_summary['absolute_bound_passed']))),
        profiles_sha256=hashlib.sha256((a.output/'profiles.npz').read_bytes()).hexdigest())
    (a.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print('passed',result['passed'],flush=True)
    if not result['passed']:raise SystemExit(1)

if __name__=='__main__':main()
