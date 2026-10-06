"""Flat scalar m6 manufactured solve: C4 factoring versus C2 endpoint constraints.

An isolated representation experiment, not a replacement binary solver.
Both use the analytic Laplacian of -2(1-t)[sqrt(t)*sqrt(1-eta^2)]^6.
A deterministic 1e-16 source perturbation probes inverse amplification.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from types import SimpleNamespace
from cartesian_modes import oracle
from numpy.polynomial import chebyshev as ch
from prolong import chebyshev_evaluation_matrix

def derivatives(n,radial):
 z=-np.cos(np.pi*(np.arange(n)+.5)/n)
 V=ch.chebvander(z,n-1);coeff=np.linalg.inv(V)
 D=np.einsum("ik,kj->ij",ch.chebvander(z,n-2),ch.chebder(coeff,axis=0),optimize=False)
 D2=np.einsum("ik,kj->ij",ch.chebvander(z,n-3),ch.chebder(coeff,m=2,axis=0),optimize=False)
 if radial:
  sigma=(1+z)/2;den=1-.8*sigma
  x=.2*sigma/den;dx=.1/den**2;ddx=.08/den**3
 else:
  th=np.tanh(2*z);x=th/np.tanh(2)
  dx=2*(1-th**2)/np.tanh(2);ddx=-8*th*(1-th**2)/np.tanh(2)
 return x,D/dx[:,None],D2/dx[:,None]**2-D*(ddx/dx**3)[:,None]

def run(n):
 nb=2*n;t,Dt,Dtt=derivatives(n,True);eta,De,Dee=derivatives(nb,False)
 tt,ee=np.meshgrid(t,eta);ss=1-ee**2;hh=4*tt/(1-tt)**2
 den=9*(hh+ss);shape=(nb,n)
 factor=lambda r:-2*(1-tt)*(tt*ss)**(r/2)
 source=factor(6)/den*(-7*(7-tt))
 rho=6*np.sqrt(tt*ss)/(1-tt)
 xyz=np.column_stack(( (3*(1+tt)/(1-tt)*ee).ravel(),rho.ravel(),np.zeros(tt.size)))
 holes=[SimpleNamespace(center=(sign*3,0,0),mass=.5) for sign in (1,-1)]
 _,_,hess,_,_=oracle(xyz,holes,3.,1.,6,False)
 independent=np.trace(hess,axis1=1,axis2=2).reshape(shape)
 oracle_error=float(np.max(abs(independent-source)/(1+abs(source))))
 assert oracle_error<1e-12
 noise=np.random.default_rng(148).standard_normal(shape)*1e-16
 # Fixed off-grid meridional points, including cylinders close to axes.
 zs=np.linspace(-.97,.97,25);ze=np.linspace(-.985,.985,27)
 ts=.2*(1+zs)/2/(1-.8*(1+zs)/2);es=np.tanh(2*ze)/np.tanh(2)
 T,E=np.meshgrid(ts,es);exact=-2*(1-T)*(T*(1-E*E))**3
 Ia=chebyshev_evaluation_matrix(n,zs);Ib=chebyshev_evaluation_matrix(nb,ze)
 rows=[]
 for r in (6,4):
  ctt=t*(1-t)**2;ct=(1-t)*((r+1)*(1-t)-2*t)
  angular=(1-eta**2)[:,None]*Dee+(-2*(r+1)*eta)[:,None]*De
  L=np.kron(np.eye(nb),ctt[:,None]*Dtt+ct[:,None]*Dt)+np.kron(angular,np.eye(n))
  diagonal=-(r+1)*(r+1-tt)+(r*r-36)*(1/hh+1/ss)
  L[np.diag_indices_from(L)]+=diagonal.ravel()
  rhs=(source*den/factor(r)).ravel();perturb=(noise*den/factor(r)).ravel()
  boundary=[]
  if r==4:
   ends=chebyshev_evaluation_matrix(nb,[-1.,1.]);start=chebyshev_evaluation_matrix(n,[-1.])[0]
   # Polar endpoint equations on two rows for every radial node. Then radial
   # endpoint equations on remaining polar rows: no duplicate corner equations.
   for j,e in ((0,ends[0]),(nb-1,ends[1])):
    for i in range(n):
     row=j*n+i;L[row]=0;L[row,i::n]=e;boundary.append(row)
   for j in range(1,nb-1):
    row=j*n;L[row]=0;L[row,j*n:(j+1)*n]=start;boundary.append(row)
   rhs[boundary]=0;perturb[boundary]=0
  assert np.isfinite(L).all() and np.isfinite(rhs).all()
  solved=np.linalg.solve(L,np.column_stack((rhs,rhs+perturb)))
  assert np.isfinite(solved).all()
  u=solved[:,0].reshape(shape);v=solved[:,1].reshape(shape)
  def sample(q):return -2*(1-T)*(T*(1-E*E))**(r/2)*np.einsum("bj,ji,ai->ba",Ib,q,Ia,optimize=False)
  clean=sample(u);noisy=sample(v)
  assert np.isfinite(clean).all() and np.isfinite(noisy).all()
  rows.append(dict(independent_cartesian_source_scaled_linf=oracle_error,formulation='factored_C4' if r==6 else 'C2_endpoint_tau',n=[n,nb],
    replaced_rows=len(boundary),manufactured_field_linf=float(np.max(abs(clean-exact))),
    noise_field_linf=float(np.max(abs(noisy-clean))),unknown_noise_linf=float(np.max(abs(v-u))),
    clean_linear_relative=float(np.linalg.norm(np.einsum("ij,j->i",L,u.ravel(),optimize=False)-rhs)/np.linalg.norm(rhs))))
 return rows

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 if a.output.exists():raise FileExistsError(a.output)
 result=dict(purpose=__doc__,driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),binary_acceptance=False,rows=[])
 for n in (12,20,32):
  result['rows']+=run(n);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['rows'][-2:]),flush=True)
