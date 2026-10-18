"""Plot the retained assessment, including a failed convergence trend."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--assessment',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
if a.output.exists():raise FileExistsError(a.output)
r=json.loads(a.assessment.read_text());rows=r['records']
assert len(rows)==3
n=[v['resolution'][0] for v in rows]
fig,axes=plt.subplots(1,2,figsize=(10,4.4),sharey=True,layout='constrained')
for ax,region,title in zip(axes,('near','bulk'),('48 near-hole points','18 bulk points')):
 for key,label,color in (('H_rms','Hamiltonian','#275f93'),('M_rms','Momentum','#bc532f')):
  ax.semilogy(n,[v['physical'][region][key] for v in rows],marker='o',label=label,color=color)
 ax.axhline(1e-6,color='#606060',linestyle='--',linewidth=1,label='Accuracy bound')
 ax.set(title=title,xlabel='Radial points (polar = 2 × radial)',xticks=n,ylim=(1e-9,2e-6))
 ax.grid(alpha=.2)
axes[0].set_ylabel(r'Physical constraint RMS [$M_0^{-2}$]')
axes[0].legend(frameon=False,fontsize=9)
fig.suptitle('One fixed QI formulation, 16 azimuthal points\nAbsolute bounds: '+('pass' if all(v['physical_bounds_passed'] for v in rows) else 'fail')+'; strict three-grid RMS trend: '+('pass' if r['decreasing_rms'] else 'fail'),fontsize=12)
fig.savefig(a.output,dpi=180);plt.close(fig)
