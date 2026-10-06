"""Plot retained evidence only; run from the repository root. No solver runs."""
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path('validation/trumpet')
OUT = Path('output/diagnostics')
OUT.mkdir(parents=True, exist_ok=True)
sources = {}
def read(name):
    p = ROOT / name
    sources[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    return json.loads(p.read_text())
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 180, 'axes.titleweight': 'bold'})
colors = ['#0072B2', '#D55E00', '#009E73', '#CC79A7']
keys = [('near','H_rms'), ('near','M_rms'), ('bulk','H_rms'), ('bulk','M_rms')]
labels = ['Near H', 'Near momentum', 'Bulk H', 'Bulk momentum']
old = read('convergence/assessment.json')['records']
tau = read('c4-experiment/axis-tau/polar-exact/sequence-v2/moderate224/result.json')
fine = read('c4-experiment/axis-tau/polar-exact/compensated-projection/gpu-v3/result.json')
spin = read('spin-focused-map/spin99_160/result.json')
boost = read('axial-factors/gamma10_512x320-v1/result.json')
fig, axs = plt.subplots(1,3,figsize=(16,5.4), layout='constrained')
for k,(reg,key) in enumerate(keys):
    axs[0].semilogy([r['resolution'][0] for r in old], [r['physical'][reg][key]/1e-6 for r in old], 'o-', color=colors[k],label=labels[k])
axs[0].set(title='Moderate spin: original sequence', xlabel=r'Radial nodes ($n_\theta=2n_r$, $n_\phi=16$)', ylabel='Physical constraint RMS / acceptance bound')
axs[0].set_xticks([192,224,256]); axs[0].legend(fontsize=9, loc='best')
for ax,records,names,title in [(axs[1],[tau,fine],['224: converged','240: stalled*'],'Moderate spin: endpoint repair'),(axs[2],[spin,boost],[r'$\chi_{seed}=0.99$'+'\n160×320×32',r'$\Gamma_{seed}=10$'+'\n512×320×8'],'Extreme binaries: retained results')]:
    for k,(reg,key) in enumerate(keys):
        ax.bar(np.arange(2)+(k-1.5)*.18,[r['physical'][reg][key]/1e-6 for r in records],width=.17,color=colors[k])
    ax.set_yscale('log'); ax.set_xticks(range(2),names); ax.set_title(title)
for ax in axs:
    ax.axhline(1,color='#333333',ls='--',lw=1.2)
    ax.grid(axis='y',alpha=.18)
fig.suptitle('Binary constraints: below the dashed line passes the RMS bound',fontsize=15)
fig.supxlabel('*240 includes compensated projection and inexact Newton; separate builds, not a controlled convergence pair.\nRMS bound = 10⁻⁶. Maxima, refinement, charges and horizon checks are additional requirements.',fontsize=10)
for ext in ['png','svg']:fig.savefig(OUT/f'trumpet_constraints.{ext}')
plt.close(fig)
h = read('convergence/horizon256/surfaces/binary.json')
fig,axs=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
xs=np.arange(len(h['records'])); ticks=['ℓ=8','ℓ=12','ℓ=16','ℓ=16\nfiner quadrature']
for i,c in enumerate(colors[:2]):
    rows=[r['holes'][i] for r in h['records']]
    chi=rows[-1]['coordinate_spin_chi']
    axs[0].semilogy(xs,[r['expansion_rms'] for r in rows],'o-',color=c,label=f'Hole {i+1}: χ ≈ {chi:.5f}')
    for key,style,label in [('mass','-','mass'),('coordinate_spin_chi','--','spin')]:
        ref=rows[-1][key]
        # Plot preceding angular-order measurements relative to independent finer quadrature.
        axs[1].semilogy(xs[:3],[abs(r[key]/ref-1) for r in rows[:3]],'o'+style,color=c,label=f'Hole {i+1} {label}')
axs[0].axhline(1e-7,color='#333333',ls='--',label='Expansion bound')
axs[1].axhline(1e-4,color='#333333',ls=':',label='Mass/spin stability bound')
axs[0].set(title='Apparent-horizon expansion',ylabel='Expansion RMS',xticks=xs,xticklabels=ticks)
axs[1].set(title='Horizon mass and spin stability',ylabel='Relative difference from finer quadrature',xticks=xs[:3],xticklabels=ticks[:3])
for ax in axs:ax.grid(axis='y',alpha=.18);ax.legend(fontsize=9)
fig.suptitle('Moderate-spin horizons work on the retained original 256×512×16 checkpoint',fontsize=13)
fig.supxlabel('Both surfaces enclose the modified interiors. These are coordinate-integral spins.\nHorizon success does not repair the binary constraint-convergence failure; no time evolution was tested here.',fontsize=10)
for ext in ['png','svg']:fig.savefig(OUT/f'trumpet_horizons.{ext}')
plt.close(fig)
(OUT/'sources.json').write_text(json.dumps({'source_sha256':sources,'generated_from_saved_results_only':True},indent=2)+'\n')
print('Saved plots and source hashes to',OUT)
