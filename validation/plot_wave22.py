"""Compact scientific summary of the one fixed Kerr packet campaign."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,required=True);p.add_argument('--figure',type=Path,required=True);a=p.parse_args()
new=json.loads((a.directory/'map12-study.json').read_text());old=json.loads((a.directory/'coarse.json').read_text());old['levels']+=json.loads((a.directory/'map2-study.json').read_text())['levels']
fig,axes=plt.subplots(1,2,figsize=(10,3.7),constrained_layout=True)
for key,label,color in [('H_rms','Hamiltonian','C0'),('M_rms','Momentum','C1')]:
 for data,style,alpha,tag in [(old,'--',.55,'focus separation 2M'),(new,'-',1.,'focus separation 12M')]:
  x=[r['n'][0] for r in data['levels']];y=[r['physical'][key] for r in data['levels']]
  axes[0].loglog(x,y,style,marker='o',alpha=alpha,color=color,label=f'{label}, {tag}')
axes[0].set(xlabel='Radial points (polar=2N, azimuthal=16)',ylabel='Unweighted physical constraint RMS',title='Fixed metric/gdot; coordinate scale matters')
axes[0].grid(alpha=.25);axes[0].legend(fontsize=7)
complex_value=lambda v:complex(*v)
raw=new['before_multipoles'];seed=new['seed_multipoles']
for key,label,color in [('metric',r'$\gamma_{22}$','C0'),('dt_metric',r'$\dot\gamma_{22}$','C1')]:
 baseline=complex_value(seed[key]);ref=complex_value(raw[key])-baseline
 values=[complex_value(r['multipoles'][key])-baseline for r in new['levels']]
 axes[1].plot([r['n'][0] for r in new['levels']],[abs(v)/abs(ref) for v in values],'-o',color=color,label=label)
axes[1].axhline(1,color='.5',ls=':',lw=1)
axes[1].set(xlabel='Radial points',ylabel='Amplitude / unsolved packet amplitude',title=r'Physical $m=2$ content at $r=10M$')
axes[1].grid(alpha=.25);axes[1].legend()
fig.suptitle(r'Liu Kerr $\chi=0.9$ + electric $\ell=m=2$ packet; constraint projection only',fontsize=11)
a.figure.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.figure,dpi=180)
