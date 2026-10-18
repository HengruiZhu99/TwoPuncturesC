"""Static plots from retained physical observers; no new solver runs."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
fig,ax=plt.subplots(figsize=(6.4,4.0),layout='constrained')
chi=np.r_[np.linspace(0,.95,200),np.linspace(.95,1,300)];q=np.sqrt(1-chi*chi)
ax.plot(chi,.5*q,label='Default quasi-isotropic')
ax.plot(chi,q,label=r'$R_0=M$ trumpet')
ax.plot(chi,(1+q)/4,label='Liu et al. wormhole')
ax.scatter([.99],[.2852669645327415],color='black',s=22,zorder=4,label=r'AthenaK Liu, $\chi=0.99$')
ax.set(xlabel=r'Isolated input spin $\chi$',ylabel=r'Coordinate horizon radius $r_H/m$',xlim=(0,1),ylim=(0,1.03),title='Coordinate charts: isolated, unboosted Kerr')
ax.text(.03,.06,r'$\chi\to1$: formal limits; exact extremality is not supported',transform=ax.transAxes,fontsize=8)
ax.legend(fontsize=8,loc='upper right');ax.grid(alpha=.18)
fig.savefig(HERE/'horizon_coordinate_radius.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(8,5.7),layout='constrained')
sets=[('Original QI map',HERE.parent/'regularity-family-control','wide','o-'),('Stable geometry',HERE/'stable-geometry','moderate','s--'),('Dual endpoint (failed)',HERE/'dual-endpoint','moderate','^-'),('Focus regularity',HERE/'focus-tau','moderate','D-')]
for label,root,prefix,style in sets:
 rows=[]
 for n in (160,192,224):
  p=root/f'{prefix}{n}'/'result.json'
  if p.exists():
   r=json.loads(p.read_text())
   if 'physical' in r:rows.append((n,r))
 for ax,(region,norm) in zip(axes.ravel(),[('near','H_rms'),('near','M_rms'),('bulk','H_rms'),('bulk','M_rms')]):
  if rows:ax.semilogy([n for n,r in rows],[r['physical'][region][norm] for n,r in rows],style,label=label)
for n in (160,192,224):
 p=HERE/'by-matched'/f'by{n}'/'result.json'
 if p.exists():
  r=json.loads(p.read_text())
  if 'physical' in r:
   for ax,region in ((axes[0,0],'near'),(axes[1,0],'bulk')):
    ax.semilogy([n],[r['physical'][region]['H_rms']],'x',color='black',ms=7,label='BY (stopping gate failed)' if n==160 else None)
for ax,(region,quantity) in zip(axes.ravel(),[('Near-hole','H'),('Near-hole','M'),('Bulk','H'),('Bulk','M')]):
 ax.axhline(1e-6,color='.5',ls=':',lw=1,label='Absolute bound' if quantity=='H' else None)
 ax.set(title=f'{region} {quantity}',xlabel=r'$n_A$ ($n_B=2n_A$, $n_\phi=16$)',ylabel=r'Physical equal-point RMS ($M_0^{-2}$)',xticks=[160,192,224]);ax.grid(alpha=.15)
axes[0,0].legend(fontsize=7)
fig.suptitle('Fixed QI free data and Cartesian observer; strict regional trend required',fontsize=10)
fig.savefig(HERE/'physical_constraints.png',dpi=180)

# One compact production comparison, on the calibrated FD6 observer.
fig,axes=plt.subplots(2,2,figsize=(8.0,5.8),layout='constrained')
for label,folder,style in [('HiSpID, first compatibility','focus-fd6','o-'),('Second compatibility','focus2-fd6','D-')]:
 rows=[]
 for n in (160,192,224):
  p=HERE/folder/f'hispid{n}'/'result.json'
  if p.exists():
   r=json.loads(p.read_text())
   if 'physical' in r:rows.append((n,r))
 for ax,(region,key) in zip(axes.ravel(),[('near','H_rms'),('near','M_rms'),('bulk','H_rms'),('bulk','M_rms')]):
  if rows:ax.semilogy([n for n,r in rows],[r['physical'][region][key] for n,r in rows],style,label=label)
for ax,region in [(axes[0,0],'near'),(axes[1,0],'bulk')]:
 rows=[(n,json.loads((HERE/'focus-fd6'/f'by{n}'/'result.json').read_text())) for n in (160,192,224)]
 ax.semilogy([n for n,r in rows],[r['physical'][region]['H_rms'] for n,r in rows],'x-',color='black',label='BY (native stopping gate failed)')
 control=json.loads((HERE/'focus-fd6/control.json').read_text())['calibration']['declared_fd6'][region]['H_rms']
 ax.axhline(control,color='.5',lw=1,ls='--',label='Exact BL observer calibration')
for ax,(region,quantity) in zip(axes.ravel(),[('Near-hole','H'),('Near-hole','M'),('Bulk','H'),('Bulk','M')]):
 ax.axhline(1e-6,color='.65',ls=':',lw=1,label='Absolute bound' if ax is axes[0,0] else None)
 ax.set(title=f'{region} {quantity}',xlabel=r'$n_A$ ($n_B=2n_A$, $n_\phi=16$)',ylabel=r'Physical equal-point RMS ($M_0^{-2}$)',xticks=[160,192,224]);ax.grid(alpha=.15)
axes[0,0].legend(fontsize=7)
fig.suptitle('Moderate QI: absolute bounds pass; strict bulk trend fails',fontsize=10)
fig.savefig(HERE/'final_constraints.png',dpi=180);plt.close(fig)
