"""Compact radial-resolution audit: collocation carrier and solved wave."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.fft import dct
from numpy.polynomial.chebyshev import chebval
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Wave22
from wave22_carrier import carrier
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,required=True);p.add_argument('--figure',type=Path,required=True);a=p.parse_args()
d=json.loads((a.directory/'physical-radial-audit.json').read_text());w=Wave22(**d['wave22'])
r=np.linspace(8,w.r0+2*np.pi*w.cycles/abs(w.omega)+4,2049)
profile_name='profiles' if (a.directory/'profiles.json').exists() else 'original-profiles'
has_profiles=(a.directory/(profile_name+'.json')).exists()
if has_profiles:
 fig=plt.figure(figsize=(8,6.5),constrained_layout=True);gs=fig.add_gridspec(2,2)
 axes=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
else:fig,axes=plt.subplots(1,2,figsize=(8,3.7),constrained_layout=True)
axes[0].plot(r,carrier(w,r)[:,0].real,'k',lw=1.5,label='Analytic carrier')
for nr,nb,sep,l,style in ((160,320,32.,.2,'--'),(576,128,16.,1.,':')):
 b=sep/2;s=(1-np.cos(np.pi*(np.arange(nr)+.5)/nr))/2;tt=l*s/(1-(1-l)*s)
 z=-np.cos(np.pi*(np.arange(nb)+.5)/nb);eta=np.tanh(2*z)/np.tanh(2.)
 nodes=b*((1+tt[:,None])/(1-tt[:,None])-eta[None,:]);co=dct(dct(carrier(w,nodes)[:,:,0][::-1,::-1],axis=0)/nr,axis=1)/nb;co[0,:]/=2;co[:,0]/=2
 C=(2*r+sep)/sep;t=(C-1)/(C+1);query=2*t/(l+(1-l)*t)-1
 from numpy.polynomial.chebyshev import chebval2d
 axes[0].plot(r,chebval2d(query,np.ones_like(r),co).real,style,label=f'{nr}x{nb}')
axes[0].set(xlabel=r'$r/M_0$',ylabel=r'$\mathrm{Re}\,q$',title='Spectral carrier representation');axes[0].legend(fontsize=9)
new=d['rows'][2:]
for key,label,color in (('value','Carrier value','C0'),('dr','Physical radial derivative','C1')):
 axes[1].loglog([v['n'][0] for v in new],[v['relative_linf'][key] for v in new],'-o',label=label,color=color)
 axes[1].plot(d['rows'][0]['n'][0],d['rows'][0]['relative_linf'][key],'x',color=color)
axes[1].axhline(1e-3,color='.4',ls=':',label='0.1% criterion');axes[1].set(xlabel='Radial collocation points',ylabel='Peak-normalized maximum error',title='Physical rays; both map coordinates');axes[1].legend(fontsize=9)
if has_profiles:
 prof=json.loads((a.directory/(profile_name+'.json')).read_text())
 with np.load(a.directory/(profile_name+'.npz')) as z:
  rr=z['radius'];fine='x'.join(map(str,prof['rows'][-1]['n']));reference=z[fine+'_ingoing'];peak=np.max(abs(reference),axis=0);angular_axes=tuple(range(1,reference.ndim))
  for row in prof['rows'][:-1]:
   key='x'.join(map(str,row['n']));coarse=z[key+'_ingoing'];error=np.max(abs(coarse-reference)/peak,axis=angular_axes)
   axes[2].semilogy(rr,np.maximum(error,1e-16),label=key)
 axes[2].axhline(1e-3,color='.4',ls=':');axes[2].set(xlabel=r'$r/M_0$',ylabel='Peak-normalized ingoing difference',title='Original solved data: ingoing mode' if profile_name=='original-profiles' else 'Solved ingoing mode');axes[2].legend(fontsize=9)
for ax in axes:ax.grid(alpha=.25)
fig.suptitle(r'Fixed $\chi_0=0.9999$ packet; radial resolution audit, no evolution',fontsize=12)
a.figure.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.figure,dpi=180)
