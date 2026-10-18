"""Read-only meridional coefficient diagnostic; no new solve or acceptance.

Compare m=0 of the factored P unknowns, not physical constraint norms.
Different families/interior windows cannot isolate a slicing-only effect.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--runs',type=Path,nargs=2,required=True)
p.add_argument('--labels',nargs=2,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--plot',type=Path,required=True)
a=p.parse_args()
if a.output.exists() or a.plot.exists():raise FileExistsError('preserve prior diagnostics')
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
def coefficients(x):
 n=x.shape[-1]
 c=np.fft.rfft(np.concatenate((x,x[...,::-1]),axis=-1),axis=-1)[...,:n]
 c=(c*np.exp(-.5j*np.pi*np.arange(n)/n)).real/n;c[...,0]*=.5
 return c

def endpoint_weights(n,side):
 # Barycentric quotient differentiated at the endpoint; trig differences
 # avoid near-endpoint subtraction. Accumulate in host extended precision.
 theta=np.arccos(np.longdouble(-1))*(np.arange(n,dtype=np.longdouble)+.5)/n
 dx=-2*np.sin(theta/2)**2 if side<0 else 2*np.cos(theta/2)**2
 value=(-1.)**np.arange(n)*np.sin(theta)/dx;value/=np.sum(value)
 derivative=value*(np.sum(value/dx)-1/dx)
 # Independent Chebyshev-polynomial endpoint values/derivatives.
 for k in (0,1,2,5):
  values=(-1.)**k*np.cos(k*theta)
  assert abs(value@values-side**k)<1e-12
  assert abs(derivative@values-(k*k*side**(k-1) if k else 0))<1e-8
 return value,derivative

rows=[]
fig,axes=plt.subplots(2,2,figsize=(10,6),layout='constrained',sharey=True)
for run,label in zip(a.runs,a.labels):
 r=json.loads((run/'result.json').read_text());artifact=run/'solve.npz'
 assert sha(artifact)==r['solve_artifact_sha256']
 assert r['library_sha256'] in (
     '7398b4b254371c39409d511ebe0808fef5e8c45b9c1e6b1af2ccb52acecf31ea',
     'c4c45eaee847a6ab17ad794c0d81321ea27d2a2cc313f84c8f0d622e7d8b4741'), 'maps must match the frozen diagnostic'
 na,nb,nphi=r['config']['n']
 assert na>=64 and nphi>=10
 with np.load(artifact) as data:all_modes=data['unknowns'].reshape(nphi,nb,na,4)
 u=all_modes[0]
 row=dict(label=label,config=r['config'],library_sha256=r['library_sha256'],
          solve_sha256=sha(artifact),result_sha256=sha(run/'result.json'),spectra={})
 for col,(axis,name) in enumerate(((1,'radial'),(0,'polar'))):
  x=np.moveaxis(u,axis,-1);n=x.shape[-1]
  c=coefficients(x)
  for line,(selection,field) in enumerate(((slice(0,1),'scalar'),(slice(1,4),'vector'))):
   amplitude=np.sqrt(np.mean(c[:,selection,:]**2,axis=(0,1)))
   scale=float(np.max(amplitude));relative=amplitude/scale
   envelope=np.array([np.max(relative[i:i+8]) for i in range(0,n,8)])
   axes[line,col].semilogy(np.arange(0,n,8)+3.5,envelope,label=label)
   axes[line,col].set(title=f'{name.capitalize()} / {field}',xlabel='Chebyshev degree',ylim=(1e-16,2))
   axes[line,col].grid(alpha=.2)
   row['spectra'][name+'_'+field]=dict(scale=scale,relative_rms_amplitude=relative.tolist(),
       last_eight_relative_rms=float(np.linalg.norm(amplitude[-8:])/np.linalg.norm(amplitude)))
 # Polynomial windows suppress one radial endpoint without changing any
 # solution. The inner endpoint includes the axis segment and both foci.
 # Exclude the final16 coefficients: degree8 multiplication can
 # alias only the top8 degrees of the retained interpolant into this DCT.
 x=np.moveaxis(u,1,-1);s=(1-np.cos(np.pi*(np.arange(na)+.5)/na))/2
 band=slice(na//2,na-16);original=coefficients(x)[...,band]
 row['radial_endpoint_suppression']={}
 for name,window in (('suppress_radial_inner_axis',s**8),('suppress_infinity',(1-s)**8)):
  c=coefficients(x*window)[...,band]
  row['radial_endpoint_suppression'][name]=dict(
      degree_band=[na//2,na-16],
      scalar_ratio=float(np.linalg.norm(c[:,:1])/np.linalg.norm(original[:,:1])),
      vector_ratio=float(np.linalg.norm(c[:,1:])/np.linalg.norm(original[:,1:])))
 # Necessary Cartesian corner compatibility for analytic low-m corrections:
 # (d_t -2 side d_eta -(m+1))P =0 at(t=0,eta=side).
 # These are representation defects, not independently measured constraints.
 # The two frozen maps are fixed at lambda=.2,kappa=2 in the retained images.
 va,da=endpoint_weights(na,-1);eta_prime=2/np.tanh(2)/np.cosh(2)**2
 row['corner_compatibility']={}
 for side in (-1,1):
  vb,db=endpoint_weights(nb,side)
  value=np.einsum('b,mbav,a->mv',vb,all_modes[:5],va)
  pt=10*np.einsum('b,mbav,a->mv',vb,all_modes[:5],da)
  pe=np.einsum('b,mbav,a->mv',db,all_modes[:5],va)/eta_prime
  defect=pt-2*side*pe-(np.arange(5)[:,None]+1)*value
  row['corner_compatibility'][str(side)]=dict(P=value.tolist(),defect=defect.tolist())
 rows.append(row)
for ax in axes[:,0]:ax.set_ylabel('Coefficient amplitude / peak')
axes[0,0].legend(frameon=False)
fig.suptitle('Retained 160 × 320 × 16 fields: factored m=0 unknowns\nDifferent seed families and interior windows; not an isolated regularity test',fontsize=11)
fig.savefig(a.plot,dpi=180);plt.close(fig)
a.output.write_text(json.dumps(dict(scope=__doc__,rows=rows,driver_sha256=sha(Path(__file__))),indent=2)+'\n')
