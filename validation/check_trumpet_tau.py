"""Verify native tau row placement/JVP and explicit checkpoint metadata."""
import argparse,json,sys,tempfile,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'python'),str(ROOT/'examples'),str(ROOT/'validation')]
from hispid import Backend
from trumpet_configs import trumpet_moderate
from checkpoint_export import write_checkpoint,read_checkpoint
from prolong import chebyshev_evaluation_matrix
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--execution',choices=['reference','kokkos'],default='reference');a=p.parse_args()
if a.output.exists():raise FileExistsError(a.output)
b=Backend(str(a.library.resolve()));assert b.parameterization().startswith('modal_P_C2tauC4_map_v5_')
c=trumpet_moderate(b,12,16);na,nb,np_=map(int,c.n);size=4*na*nb*np_;k=np.arange(size)
u=1e-3*np.sin(.13*k);d=1e-3*np.cos(.27*k)
phi=2*np.pi*np.arange(np_)/np_;basis=np.zeros((np_,np_))
for slot in range(np_):
 m=slot if slot<=np_//2 else slot-np_//2
 basis[:,slot]=np.sqrt((1 if m in (0,np_//2) else 2)/np_)*(np.cos(m*phi) if slot<=np_//2 else np.sin(m*phi))
with b.create(c,execution=a.execution,geometry='host') as s:
 s.set_unknowns(u);actual=s.residual(u);eq=s.equation_samples();raw=eq['physical_equivalent'].copy()
 raw[:,0]*=-eq['psi']**5/8;raw[:,1:]*=eq['psi'][:,None]**10
 weight=(np.sin(np.pi*(np.arange(nb)+.5)/nb)[:,None]*np.sin(np.pi*(np.arange(na)+.5)/na)[None,:])**3
 raw=raw.reshape(np_,nb,na,4)*weight[None,:,:,None]
 expected=np.einsum('km,kjiv->mjiv',basis,raw,optimize=False)
 values=u.reshape(np_,nb,na,4);ends=chebyshev_evaluation_matrix(nb,[-1.,1.]);start=chebyshev_evaluation_matrix(na,[-1.])[0]
 for slot in range(np_):
  m=slot if slot<=np_//2 else slot-np_//2
  if m<5:continue
  expected[slot,0]=np.einsum('j,jiv->iv',ends[0],values[slot],optimize=False)
  expected[slot,-1]=np.einsum('j,jiv->iv',ends[1],values[slot],optimize=False)
  expected[slot,1:-1,0]=np.einsum('i,jiv->jv',start,values[slot,1:-1],optimize=False)
 actual_modes=np.einsum('km,kjiv->mjiv',basis,actual.reshape(np_,nb,na,4),optimize=False)
 row_error=float(np.max(abs(expected-actual_modes)))
 j=s.jvp(u,d);h=.005;fd=(s.residual(u+h*d)-s.residual(u-h*d))/(2*h)
 jvp_error=float(np.linalg.norm(fd-j)/np.linalg.norm(j))
with tempfile.TemporaryDirectory() as tmp:
 path=Path(tmp)/'tau.checkpoint';write_checkpoint(path,c,u,b.loaded_sha256,'diagnostic',b.parameterization());_,v,meta=read_checkpoint(path)
 assert np.array_equal(v,u) and meta['parameterization']==b.parameterization()
result=dict(row_error_linf=row_error,jvp_relative_l2=jvp_error,checkpoint_roundtrip=True,library_sha256=b.loaded_sha256,driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),execution=a.execution,binary_acceptance=False)
assert row_error<1e-12 and jvp_error<1e-8
result['passed']=True;a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
