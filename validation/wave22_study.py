"""One fixed Liu Kerr + SpEC electric packet, three production-CUDA levels.

Raw physical constraints use independent FD6 Cartesian stencils at fixed
observers. Norms are unweighted point RMS, not the native scaled residual.
This is initial-data validation, not a superradiant evolution experiment.
"""
import argparse
import ctypes as C
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Backend,Wave22
from kerr_wave_case import config,gauge
from physical import constraints,norms
from parallel_observer import ParallelObserver
from checkpoint_export import read_checkpoint


def multipoles(solution,radius=10.,threads=16,chi=.9):
    mu,wt=np.polynomial.legendre.leggauss(16);phi=2*np.pi*np.arange(24)/24
    mu,phi=np.meshgrid(mu,phi,indexing='ij');sn=np.sqrt(1-mu*mu);cp=np.cos(phi);sp=np.sin(phi)
    n=np.stack((sn*cp,sn*sp,mu),axis=-1).reshape(-1,3);x=radius*n
    et=np.stack((mu*cp,mu*sp,-sn),axis=-1).reshape(-1,3)
    ep=np.stack((-sp,cp,np.zeros_like(mu)),axis=-1).reshape(-1,3)
    with ParallelObserver(solution,threads) as sample:q=sample.with_derivatives(x)
    g=q['gamma'].reshape(-1,3,3);K=q['Kij'].reshape(-1,3,3)
    alpha,beta=gauge(x,chi);db=[];eps=1e-4
    for d in range(3):db.append((gauge(x+eps*np.eye(3)[d],chi)[1]-gauge(x-eps*np.eye(3)[d],chi)[1])/(2*eps))
    db=np.array(db);dg=q['dgamma']
    dtg=np.einsum('nk,nkij->nij',beta,dg)+np.einsum('nkj,ink->nij',g,db)+np.einsum('nik,jnk->nij',g,db)-2*alpha[:,None,None]*K
    Y=np.sqrt(5/(64*np.pi))*(1+mu.ravel())**2*np.exp(2j*phi.ravel())
    weight=np.repeat(wt,24)*(2*np.pi/24)
    out={}
    for name,T in (('metric',g),('dt_metric',dtg),('dr_metric',np.einsum('nk,nkij->nij',n,dg))):
        polarization=.5*(np.einsum('ni,nij,nj->n',et,T,et)-np.einsum('ni,nij,nj->n',ep,T,ep))-1j*np.einsum('ni,nij,nj->n',et,T,ep)
        a=np.sum(weight*np.conj(Y)*polarization);out[name]=[float(a.real),float(a.imag)]
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--levels',type=int,nargs='+',default=[32,48,64])
    p.add_argument('--focus-separation',type=float,default=2.);p.add_argument('--execution',default='kokkos');p.add_argument('--threads',type=int,default=16)
    p.add_argument('--initial-checkpoint',type=Path,help='resume the first matching grid from an explicitly retained checkpoint')
    p.add_argument('--polish-linear-rtol',type=float,default=0.,help='retry a near-converged failed inner solve with inexact forcing, retaining outer tolerance')
    p.add_argument('--max-krylov',type=int,default=1200)
    p.add_argument('--max-newton',type=int,default=12)
    p.add_argument('--krylov-restart',type=int,default=100)
    p.add_argument('--krylov',choices=['gmres','lgmres','bicgstab'],default='gmres')
    p.add_argument('--linear-rtol',type=float,default=0.,help='inexact Newton linear tolerance; outer config tolerance stays fixed')
    p.add_argument('--polar-points',type=int,help='hold polar resolution fixed for a radial-only sequence')
    p.add_argument('--memory-limit-mib',type=int,default=32768)
    p.add_argument('--chi',type=float,default=.9);p.add_argument('--amplitude',type=float,default=.001)
    f=p.add_mutually_exclusive_group();f.add_argument('--omega',type=float);f.add_argument('--frequency-ratio',type=float,help='carrier / (m Omega_H), m=2, mass=1')
    p.add_argument('--cycles',type=float,default=1.);p.add_argument('--r0',type=float,default=6.)
    p.add_argument('--multipole-radius',type=float,default=10.)
    p.add_argument('--extended-observers',action='store_true',help='include near-horizon and packet-wide radii')
    a=p.parse_args()
    if not 0<a.chi<1:raise ValueError('positive subextremal spin required')
    if len(a.levels)<1 or any(n<8 for n in a.levels) or any(j<=i for i,j in zip(a.levels,a.levels[1:])):raise ValueError('increasing levels >=8 required')
    threshold=a.chi/(1+np.sqrt(1-a.chi*a.chi))
    omega=a.frequency_ratio*threshold if a.frequency_ratio is not None else (.5 if a.omega is None else a.omega)
    a.output.mkdir(parents=True,exist_ok=False);b=Backend(str(a.library.resolve()))
    if not np.isfinite(a.focus_separation) or a.focus_separation<=0:raise ValueError('positive focus separation required')
    def case_config(n=8):
        c=config(b,n,a.chi);c.hole[1].center[:]=(0,0,-a.focus_separation)
        if a.polar_points is not None:c.n[1]=a.polar_points
        c.memory_limit_mib=a.memory_limit_mib;c.max_krylov=a.max_krylov;c.max_newton=a.max_newton;c.krylov_restart=a.krylov_restart
        return c
    w=Wave22(amplitude=a.amplitude,omega=omega,cycles=a.cycles,r0=a.r0,cutoff_min=2.,cutoff_max=4.)
    length=2*np.pi*w.cycles/abs(w.omega)
    horizon_radius=(1+np.sqrt(1-a.chi*a.chi))/4
    dirs=np.array([[.3,.4,.5],[-.7,.2,.3],[.2,-.8,.4]]);dirs/=np.linalg.norm(dirs,axis=1)[:,None]
    radii=np.array([.5,.8,1.5,3.,5.,6.5,9.,13.,17.,20.])
    if a.extended_observers:radii=np.unique(np.concatenate((radii,[1.1*horizon_radius,1.5*horizon_radius,w.r0,w.r0+.25*length,w.r0+.5*length,w.r0+.75*length,w.r0+length,w.r0+length+3.])))
    xyz=(radii[:,None,None]*dirs[None,:,:]).reshape(-1,3)
    steps=np.minimum(.005,.006*np.linalg.norm(xyz,axis=1))
    report=dict(formulation='fixed physical gamma/gdot perturbation, stationary Liu gauge, coupled CTT',
        chi=a.chi,mass=1.,seed_horizon_radius=horizon_radius,seed_superradiant_threshold=threshold,carrier_threshold_ratio=omega/threshold,multipole_radius=a.multipole_radius,inactive_focus_separation=a.focus_separation,wave22=w.as_dict(),execution=a.execution,geometry='host-long-double',
        library_sha256=b.library_sha256(),parameterization=b.parameterization(),parameterization_maps=b.parameterization_maps(),krylov=a.krylov,max_krylov=a.max_krylov,krylov_restart=a.krylov_restart,max_newton=a.max_newton,linear_rtol=a.linear_rtol,polish_linear_rtol=a.polish_linear_rtol,
        source_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in (Path(__file__),ROOT/'src/HiSpID_free_data.hpp',ROOT/'validation/physical.py',ROOT/'validation/kerr_wave_case.py')},
        observer=dict(xyz=xyz.tolist(),steps=steps.tolist(),order=6,norm='unweighted fixed-point RMS'),levels=[])
    def save():(a.output/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    with b.create_sampler(case_config()) as seed,b.create_sampler(case_config(),wave=w) as trial:
        with ParallelObserver(seed,a.threads) as sample:report['exact_seed_floor']=norms(constraints(sample,xyz,steps,6))
        with ParallelObserver(trial,a.threads) as sample:report['before_solve']=norms(constraints(sample,xyz,steps,6))
        report['seed_multipoles']=multipoles(seed,a.multipole_radius,a.threads,a.chi);report['before_multipoles']=multipoles(trial,a.multipole_radius,a.threads,a.chi)
    save();print('before',report['before_solve'],flush=True)
    initial_values=None;initial_metadata=None
    if a.initial_checkpoint:
        cc,initial_values,initial_metadata=read_checkpoint(a.initial_checkpoint)
        def signature(c):
            out={}
            for name,_ in c._fields_:
                if name in ('max_newton','max_krylov','krylov_restart','memory_limit_mib'):continue
                v=getattr(c,name)
                out[name]=[[h.mass,*h.center,*h.spin,*h.velocity] for h in v] if name=='hole' else (list(v) if isinstance(v,C.Array) else v)
            return out
        if signature(cc)!=signature(case_config(a.levels[0])) or cc.wave22.as_dict()!=w.as_dict() or initial_metadata['parameterization']!=b.parameterization() or cc.seed_family!='liu':raise ValueError('initial checkpoint differs from the fixed first-grid case')
    previous=None
    for index,n in enumerate(a.levels):
        start=time.monotonic();c=case_config(n)
        with b.create(c,execution=a.execution,wave=w) as s:
            resumed=initial_values is not None and index==0
            if resumed:s.set_unknowns(initial_values)
            initial=float(np.max(abs(s.residual(s.unknowns()))));print("setup",n,s.setup_statistics(),flush=True)
            first_forcing=a.polish_linear_rtol if resumed and a.polish_linear_rtol else a.linear_rtol
            diag=s.solve(krylov=a.krylov,linear_rtol=first_forcing)
            attempts=[dict(diagnostics=diag,linear_rtol=first_forcing,linear_history=s.linear_history())]
            if a.polish_linear_rtol and not diag['converged'] and max(diag['scaled_linf'])<.1*initial and s.linear_history() and s.linear_history()[-1][2]<.5:
                print('decreased residual; retaining outer tolerance and continuing inexact Newton',n,diag,flush=True)
                diag=s.solve(krylov=a.krylov,linear_rtol=a.polish_linear_rtol)
                attempts.append(dict(diagnostics=diag,linear_rtol=a.polish_linear_rtol,linear_history=s.linear_history()))
            print('solved',n,diag,flush=True);u=s.unknowns();fields=s.sample(xyz)
            row=dict(n=list(c.n),initial_scaled_linf=initial,initial_checkpoint=initial_metadata if resumed else None,diagnostics=diag,solve_attempts=attempts,total_solve_seconds=sum(v['diagnostics']['seconds'] for v in attempts),setup=s.setup_statistics(),
                     work=s.work_statistics(),linear_history=s.linear_history())
            with ParallelObserver(s,a.threads) as sample:row['physical']=norms(constraints(sample,xyz,steps,6))
            row['multipoles']=multipoles(s,a.multipole_radius,a.threads,a.chi)
            if previous is not None:row['previous_solution_difference']={k:float(np.max(abs(fields[k]-previous[k]))) for k in ('gamma','Kij')}
            np.savez_compressed(a.output/f'n{n}.npz',unknowns=u,xyz=xyz,**{k:fields[k] for k in ('gamma','Kij','correction')})
            row['wall_seconds']=time.monotonic()-start
            row['checkpoint']=s.write_checkpoint(a.output/f'n{n}.hispid','preliminary' if diag['converged'] else 'diagnostic')
            previous=fields
            report['levels'].append(row);save();print('level',n,row,flush=True)
            if not diag['converged']:raise RuntimeError('nonconverged solve retained for diagnosis')
    if len(report['levels'])==3:
        records=report['levels'];report['assessment']={
            'all_native_converged':all(r['diagnostics']['converged'] for r in records),
            'hamiltonian_rms_decreases':all(records[i+1]['physical']['H_rms']<records[i]['physical']['H_rms'] for i in (0,1)),
            'momentum_rms_decreases':all(records[i+1]['physical']['M_rms']<records[i]['physical']['M_rms'] for i in (0,1)),
            'physical_solution_differences_decrease':{k:records[2]['previous_solution_difference'][k]<records[1]['previous_solution_difference'][k] for k in ('gamma','Kij')},
            'note':'Finite fixed observer set, no global bound, evolution or measured superradiant amplification claimed.'}
        def characteristic(v):
            z={k:complex(*v[k])-complex(*report['seed_multipoles'][k]) for k in ('metric','dt_metric','dr_metric')}
            dt=a.multipole_radius*z['dt_metric'];dr=z['metric']+a.multipole_radius*z['dr_metric']
            incoming=(dt+dr)/2;outgoing=(dt-dr)/2
            return dict(ingoing_amplitude=abs(incoming),outgoing_amplitude=abs(outgoing),outgoing_to_ingoing=abs(outgoing)/abs(incoming))
        report['asymptotic_metric_characteristics']=dict(radius=a.multipole_radius,
            definition='H=r delta(gamma22); ingoing=(dtH+drH)/2, outgoing=(dtH-drH)/2; coordinate diagnostic, not curvature energy flux',
            before=characteristic(report['before_multipoles']),after=[characteristic(r['multipoles']) for r in records])
        save()

if __name__=='__main__':main()
