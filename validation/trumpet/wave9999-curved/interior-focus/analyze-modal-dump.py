"""Spatial/spectral localization of the exact retained base-inverse response."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.fft import dctn
p=argparse.ArgumentParser();p.add_argument('prefix',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
na,nb,np_=384,128,16
shape=(np_,nb,na,4)
paths={k:Path(str(a.prefix)+'.'+k+'.f64') for k in ('rhs','correction','scale')}
v={k:np.fromfile(path,dtype='<f8').reshape(shape) for k,path in paths.items()}
assert all(np.isfinite(x).all() for x in v.values())
# Verify the axis labels and degree indexing with one analytic product.
witness=np.cos(3*np.pi*(np.arange(nb)+.5)/nb)[:,None]*np.cos(7*np.pi*(np.arange(na)+.5)/na)[None,:]
w=dctn(witness,norm='ortho');assert np.unravel_index(np.argmax(abs(w)),w.shape)==(3,7)
w[3,7]=0;assert np.linalg.norm(w)<1e-10
phi=2*np.pi*np.arange(np_)/np_
F=np.array([np.sqrt((1 if m in (0,np_//2) else 2)/np_)*(np.cos(m*phi) if slot<=np_//2 else np.sin(m*phi)) for slot in range(np_) for m in [slot if slot<=np_//2 else slot-np_//2]])
modal_rhs=np.einsum('kl,lijv->kijv',F.astype(np.longdouble),v['rhs'].astype(np.longdouble)).astype(float)
weight=(np.sin(np.pi*(np.arange(nb)+.5)/nb)[:,None]*np.sin(np.pi*(np.arange(na)+.5)/na)[None,:])**3
physical_rows=modal_rhs/weight[None,:,:,None]
# Focus rows are regularity equations, not weighted PDE samples. Exclude them.
for slot in range(np_):
 mode=slot if slot<=np_//2 else slot-np_//2
 if mode<=4: physical_rows[slot,[0,nb-1],0,:]=0
 else:
  physical_rows[slot,:,0,:]=0;physical_rows[slot,[0,nb-1],:,:]=0
load=modal_rhs/v['scale']
sig=.5*(1-np.cos(np.pi*(np.arange(na)+.5)/na));t=16*sig/(1+15*sig);C=(1+t)/(1-t)
zeta=-np.cos(np.pi*(np.arange(nb)+.5)/nb);eta=np.tanh(.3*zeta)/np.tanh(.3)
r=.5*(C[None,:]-eta[:,None]);rho=.5*np.sqrt((C[None,:]**2-1)*(1-eta[:,None]**2))
inner=r<.1;axis=np.zeros((nb,na),bool);axis[:4]=True;axis[-4:]=True;radial=np.zeros_like(axis);radial[:,:4]=True;outer=np.zeros_like(axis);outer[:,-4:]=True
rows=[]
for mode in range(5):
 slots=[mode]+([np_//2+mode] if mode else [])
 entry={'mode':mode}
 for name,array in [('correction',v['correction']),('rhs_modal',modal_rhs),('pde_rhs_weight_removed',physical_rows),('inverse_load',load)]:
  field=array[slots,:,:,1:];spec=dctn(field,type=2,axes=(1,2),norm='ortho');energy=np.sum(field*field);en=spec*spec
  maxindex=np.unravel_index(np.argmax(abs(field)),field.shape);_,j,i,comp=maxindex
  entry[name]=dict(l2=float(np.sqrt(energy)),maximum=float(np.max(abs(field))),maximum_i=int(i),maximum_j=int(j),maximum_component=int(comp+1),maximum_radius=float(r[j,i]),maximum_rho=float(rho[j,i]),polar_upper_half_fraction=float(np.sum(en[:,nb//2:,:,:])/energy),radial_upper_half_fraction=float(np.sum(en[:,:,na//2:,:])/energy),polar_last8_fraction=float(np.sum(en[:,-8:,:,:])/energy),radial_last8_fraction=float(np.sum(en[:,:,-8:,:])/energy),near_polar_edge_fraction=float(np.sum(field[:,axis,:]**2)/energy),near_inner_radial_edge_fraction=float(np.sum(field[:,radial,:]**2)/energy),near_outer_radial_edge_fraction=float(np.sum(field[:,outer,:]**2)/energy),interior_fraction=float(np.sum(field[:,inner,:]**2)/energy),parseval_relative=float(abs(np.sum(en)-energy)/energy))
 rows.append(entry)
out=dict(scope='Coordinate angular-edge masks do not imply physical proximity to the axis. PDE RHS with native weight removed excludes boundary equations. Unweighted modal P and modal/projected RHS spectra, not physical constraint norms. Original momentum RHS normalized to one; exact same frozen five-update checkpoint.',grid=[na,nb,np_],files={k:dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for k,path in paths.items()},rows=rows)
a.output.write_text(json.dumps(out,indent=2)+'\n')
for row in rows:print(row)
