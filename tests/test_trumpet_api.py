"""Compact seed-family integration controls; one native production candidate.

Run: python3 tests/test_trumpet_api.py --library /absolute/libHiSpID.so
"""
import argparse
import ctypes as C
import json
import sys
import tempfile
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
from hispid import Backend,Hole,Config
from checkpoint_export import read_checkpoint,write_checkpoint

def run(path):
    b=Backend(str(path.resolve()));c=b.config()
    c.n[:]=[8,8,4];c.far_radius=0;c.inner_flatten=0;c.inner_max[:]=[0,0]
    c.hole[0]=Hole(1,center=(0,0,0),spin=(.2,-.3,.7),velocity=(.25,.1,-.2))
    c.hole[1].mass=0;c.seed_family='trumpet_r0_m'
    x=np.array([[1.2,.7,-.4],[2,-1,.8],[4,.3,2.]])
    seed=b.seed(c.hole[0],x,seed_family=c.seed_family)
    with b.create_sampler(c) as s:
        assert b.lib.HiSpID_seed_family(s.context)==1
        sampled=s.sample(x)
        err=max(float(np.max(np.abs(sampled[k]-seed[k])/(1+np.abs(seed[k])))) for k in ['gamma','Kij'])
        assert err<1e-12
    # A generic rigid rotation and mass scaling are distinct convention checks.
    q,_=np.linalg.qr(np.array([[.3,.5,.2],[.7,-.2,.4],[.1,.8,-.9]]))
    if np.linalg.det(q)<0:q[:,0]*=-1
    shift=np.array([.2,-.7,1.3]);mass=2.3
    h=Hole(mass,center=shift,spin=mass**2*(q@np.array(c.hole[0].spin)),velocity=q@np.array(c.hole[0].velocity))
    transformed=b.seed(h,mass*x@q.T+shift,seed_family=c.seed_family)
    covariance=[]
    for key,factor in [('gamma',1),('Kij',1/mass)]:
        expected=factor*np.einsum('ia,nab,jb->nij',q,seed[key].reshape(-1,3,3),q)
        covariance.append(float(np.max(np.abs(expected-transformed[key].reshape(-1,3,3))/(1+np.abs(expected)))))
    assert max(covariance)<1e-11
    # Native JVP at a nonzero state, two step sizes; no nonlinear solve.
    with b.create(c) as s:
        z=np.zeros(s.size);seed_res=float(np.max(np.abs(s.residual(z))))
        assert seed_res<1e-10
        rng=np.random.default_rng(721);base=rng.normal(0,1e-6,s.size);direction=rng.normal(0,1e-3,s.size)
        exact=s.jvp(base,direction);jvp=[]
        for step in [1e-3,5e-4]:
            fd=(s.residual(base+step*direction)-s.residual(base-step*direction))/(2*step)
            jvp.append(float(np.linalg.norm(fd-exact)/np.linalg.norm(exact)))
        assert max(jvp)<1e-7
    with tempfile.TemporaryDirectory() as tmp:
        p=Path(tmp)/'trumpet.checkpoint';unknowns=np.zeros(4*np.prod(c.n))
        write_checkpoint(p,c,unknowns,b.loaded_sha256,'analytic_seed')
        decoded,values,meta=read_checkpoint(p)
        assert meta['format_version']==2 and decoded.seed_family==c.seed_family
        with b.create_sampler(decoded) as s:assert b.lib.HiSpID_seed_family(s.context)==1
        original=p.read_text();rejected=0
        for altered in [original.replace('trumpet_r0_m','other'),original.replace('HISPID_CHECKPOINT 2','HISPID_CHECKPOINT 1'),original.replace('seed_family trumpet_r0_m\n','')]:
            p.write_text(altered)
            try:read_checkpoint(p)
            except ValueError:rejected+=1
        assert rejected==3
        qi=Config.from_buffer_copy(c)
        write_checkpoint(p,qi,unknowns,b.loaded_sha256,'analytic_seed')
        qc,qv,qm=read_checkpoint(p)
        assert qm['format_version']==1 and qc.seed_family=='qi'
        with b.create_sampler(qc) as s:assert b.lib.HiSpID_seed_family(s.context)==0
    assert not b.lib.HiSpID_create_with_seed_family(C.byref(c),17,0,0,1)
    return dict(sampler_scaled_difference=err,covariance_scaled_difference=covariance,
                seed_weighted_residual=seed_res,jvp_relative_l2=jvp,
                checkpoint_rejection_cases=3,checkpoint_versions=[1,2],
                native_library_sha256=b.loaded_sha256,passed=True,binary_acceptance=False)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--library',type=Path,required=True)
    ap.add_argument('--output',type=Path);args=ap.parse_args()
    if args.output and args.output.exists():raise FileExistsError(args.output)
    result=run(args.library);text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    print(text)
