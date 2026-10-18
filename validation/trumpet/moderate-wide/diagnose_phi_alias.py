"""Double azimuthal sampling of a frozen field, without solving or acceptance."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Backend
from checkpoint_export import read_checkpoint
from physical import constraints
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--checkpoint',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
if a.output.exists():raise FileExistsError(a.output)
c,u,meta=read_checkpoint(a.checkpoint);b=Backend(str(a.library.resolve()));assert b.loaded_sha256==meta['source_library_sha256']
na,nb,oldphi=map(int,c.n);assert oldphi==16;np_=2*oldphi;half=oldphi//2
old=u.reshape(oldphi,nb,na,4);values=np.zeros((np_,nb,na,4))
for m in range(half+1):
 scale=np.sqrt((1 if m in (0,half) else 2)/oldphi)/np.sqrt((1 if m==0 else 2)/np_)
 values[m]=scale*old[m]
 if m and m!=half:values[np_//2+m]=scale*old[half+m]
c.n[2]=np_;c.krylov_restart=8;c.memory_limit_mib=65536
started=time.monotonic();a.output.parent.mkdir(parents=True,exist_ok=True)
out=dict(scope=__doc__,checkpoint=meta,library_sha256=b.loaded_sha256,shape=list(c.n),binary_acceptance=False)
def save():a.output.write_text(json.dumps(out,indent=2)+'\n')
save()
with b.create(c,execution='kokkos',geometry='host') as s:
 out['setup_seconds']=time.monotonic()-started;save();s.set_unknowns(values.ravel());eq=s.equation_samples()
 x=eq['xyz'].reshape(np_,nb,na,3);hm=eq['physical_equivalent'].reshape(np_,nb,na,4)
 rad=np.linalg.norm(x[0],axis=2);rho=np.linalg.norm(x[0,:,:,1:],axis=2)
 mask=(rad>=4)&(rad<=16)&(rho>1)&np.all(eq['attenuation'].reshape(np_,nb,na)==1,axis=0)
 mask[[0,-1]]=False;mask[:,[0,-1]]=False;assert mask.sum()>50
 out['meridional_points']=int(mask.sum())
 for label,sl in [('old_phases',slice(None,None,2)),('interleaved_phases',slice(1,None,2))]:
  v=hm[sl][:,mask,:];out[label]=dict(component_rms=np.sqrt(np.mean(v*v,axis=(0,1))).tolist(),component_max=np.max(abs(v),axis=(0,1)).tolist())
 f=np.fft.rfft(hm[:,mask,:],axis=0)/np_;factor=np.full(len(f),2.);factor[[0,-1]]=1
 energy=factor[:,None,None]*abs(f)**2
 out['azimuthal_component_energy_fraction']=(energy.sum(axis=1)/energy.sum(axis=(0,1))).tolist()
 out['cosine_component_rms']=np.sqrt(np.mean((f.real*factor[:,None,None])**2,axis=1)).tolist()
 out['sine_component_rms']=np.sqrt(np.mean((f.imag*factor[:,None,None])**2,axis=1)).tolist();save()
 # Confirm selected interleaved samples with the independent physical observer.
 score=np.linalg.norm(hm[1::2][:,mask,1:],axis=2);ids=np.argsort(score.ravel())[-4:]
 odd,mer=np.unravel_index(ids,score.shape);sel=np.flatnonzero(mask.ravel())[mer]
 xx=x.reshape(np_,nb*na,3)[2*odd+1,sel];expected=hm.reshape(np_,nb*na,4)[2*odd+1,sel]
 dist=np.min(np.linalg.norm(xx[:,None,:]-np.array([h.center for h in c.hole])[None,:,:],axis=2),axis=1)
 obs=constraints(s.sample,xx,np.minimum(.002,.002*dist))
 measured=np.column_stack([obs['H'],obs['M']]);out['independent_witness']=dict(xyz=xx.tolist(),native=expected.tolist(),physical_fd=measured.tolist(),difference_linf=float(np.max(abs(expected-measured))))
 out['seconds']=time.monotonic()-started;save()
print(json.dumps({k:out[k] for k in ['old_phases','interleaved_phases','independent_witness','seconds']},indent=2),flush=True)
