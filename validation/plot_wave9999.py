"""One compact diagnostic figure for the fixed near-threshold Kerr packet."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--study',type=Path,required=True);p.add_argument('--figure',type=Path,required=True)
a=p.parse_args();d=json.loads(a.study.read_text());rows=d['levels'];x=[r['n'][0] for r in rows]
fig,axes=plt.subplots(1,3,figsize=(11,3.5),constrained_layout=True)
for key,label,color in [('H_rms','Hamiltonian','C0'),('M_rms','Momentum','C1')]:
 axes[0].semilogy(x,[r['physical'][key] for r in rows],'-o',label=label,color=color)
 axes[0].axhline(d['exact_seed_floor'][key],ls=':',color=color,lw=1)
axes[0].set(ylabel='Physical constraint RMS',title='Fixed exterior observers\nDotted: unperturbed observer floor')
axes[0].legend(fontsize=8)
for key,label in [('gamma',r'$\max|\Delta\gamma_{ij}|$'),('Kij',r'$\max|\Delta K_{ij}|$')]:
 axes[1].semilogy(x[1:],[r['previous_solution_difference'][key] for r in rows[1:]],'-o',label=label)
axes[1].set(ylabel='Difference from preceding level',title='Sampled physical solution');axes[1].legend(fontsize=8)
char=d['asymptotic_metric_characteristics'];ratios=[v['outgoing_to_ingoing'] for v in char['after']]
axes[2].plot(x,ratios,'-o',color='C2')
axes[2].set(ylabel=r'$|C_{\mathrm{out}}|/|C_{\mathrm{in}}|$',title=f"Metric characteristics at r={char['radius']:.2f}M\nAfter constraint projection")
axes[2].ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
for ax in axes:ax.set_xlabel('Radial points (polar=2N; azimuthal=16)');ax.grid(alpha=.25);ax.set_xticks(x)
fig.suptitle(r'Liu Kerr $\chi_{\rm seed}=0.9999$; ingoing prograde packet, $\omega=0.99(2\Omega_H)$',fontsize=11)
a.figure.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.figure,dpi=180)
