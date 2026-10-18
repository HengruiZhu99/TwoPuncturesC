"""Liu--Etienne--Shapiro chart: independent radial pullback and vacuum checks.

This validates exact seeds, not a superposed binary solve or NR evolution.
Reference: arXiv:1001.4077, Eqs.12--15 (PRD80,121503).
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Backend,Hole
from physical import constraints,norms

def check(library,output,execution='reference'):
    b=Backend(str(library.resolve()));rows=[];worst=0
    # Include the throat itself: signed lapse/sqrt(Delta) cancellation must
    # remain finite, as well as both sides of the bridge and a generic axis.
    for chi in (0.,.39,.99):
        rp=1+np.sqrt(1-chi*chi);rm=chi*chi/rp;c=rp/4
        directions=np.array([[.3,.4,.5],[.7,-.2,.4],[-.5,.1,.3],[0,0,1.]])
        directions/=np.linalg.norm(directions,axis=1)[:,None]
        radius=np.array([.7*c,c,1.7*c,2.3]);x=directions*radius[:,None]
        rb=radius+2*c+c*c/radius;root=np.sqrt((rb-rm)/radius)
        signed=(radius-c)*root;rbprime=1-c*c/radius**2
        rootprime=((rbprime*radius-(rb-rm))/radius**2)/(2*root)
        q=(rb-1+signed)/2;qp=(rbprime+root+(radius-c)*rootprime)/2
        jac=q[:,None,None]/radius[:,None,None]*np.eye(3)+np.einsum('n,ni,nj->nij',qp-q/radius,directions,directions)
        h=Hole(1,spin=(0,0,chi));ref=b.seed(h,directions*q[:,None],choice=0)
        got=b.seed(h,x,choice=0,seed_family='liu',execution=execution)
        differences={}
        for field in ('gamma','Kij'):
            expected=np.swapaxes(jac,1,2)@ref[field].reshape(-1,3,3)@jac
            actual=got[field].reshape(-1,3,3)
            differences[field]=float(np.max(abs(actual-expected)/(1+abs(expected))))
        worst=max(worst,*differences.values());assert max(differences.values())<2e-12,differences
        # Independent Cartesian finite differences, at fixed physical points.
        observer=np.array([[.42,.3,.2],[.1,-.25,.05],[1.4,-.8,.6],[4.,2.,-3.]])
        for velocity in ((0.,0.,0.),(.2,-.1,.05)):
            h.velocity[:]=velocity
            sample=lambda y:b.seed(h,y,choice=0,seed_family='liu',execution=execution)
            checks=[dict(step=step,norms=norms(constraints(sample,observer,step))) for step in (.002,.001)]
            assert checks[-1]['norms']['H_max']<1e-7 and checks[-1]['norms']['M_max']<1e-7,checks
            rows.append(dict(chi=chi,velocity=velocity,horizon_radius=c,pullback=differences,vacuum=checks))
        # Minimal API/sampler check: one exact seed + zero correction.
        cfg=b.config();cfg.conformal_choice=0;cfg.seed_family='liu';cfg.hole[0]=h;cfg.hole[1].mass=0;cfg.inner_min[:]=[0,0];cfg.inner_max[:]=[0,0];cfg.omega[:]=[0,0];cfg.far_radius=0
        with b.create_sampler(cfg) as s:
            z=s.sample(observer);exact=sample(observer)
            for field in ('gamma','Kij','psi'):
                assert np.max(abs(z[field]-exact[field])/(1+abs(exact[field])))<2e-12
    result=dict(passed=True,scope='exact_seeds_only',library_sha256=b.library_sha256(),execution=execution,
                worst_radial_pullback_scaled_linf=worst,cases=rows,
                extremal_limit_horizon_radius_over_M=.25,
                reference='https://arxiv.org/html/1001.4077#S2.SS1')
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='cases'}))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--execution',choices=('reference','kokkos'),default='reference');a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    check(a.library,a.output,a.execution)
