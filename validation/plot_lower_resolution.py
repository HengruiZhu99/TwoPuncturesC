"""Concise solution/constraint diagnostics from saved results; no solves."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--latex-figures', type=Path)
    a = p.parse_args()
    root = a.results
    load = lambda path:json.loads((root/path).read_text())
    moderate = load(Path('moderate/result.json'))
    prefix = load(Path('moderate-prefix/result.json'))
    single = load(Path('single-hole-final/result.json'))
    native = load(Path('single-native.json'))
    pullback = load(Path('liu-pullback-9999.json'))
    paper = load(Path('paper-binary-reference-norm-bound.json'))
    rows = sorted(prefix['records']+moderate['records'], key=lambda r:r['n'][0])
    checks = dict(single_hole_analytic_liu=single['liu_passed'],
                  independent_liu_pullback=pullback['passed'], native_single_hole=native['passed'],
                  same_equation_controls=[], physical_free_data_fixed=True,
                  high_resolution_bulk_issue_resolved=False,
                  paper_reference_native_gate=load(Path('paper-binary/result.json'))['records'][-2]['diagnostics']['converged'])
    excluded = {'n', 'tolerance', 'max_newton', 'max_krylov', 'krylov_restart', 'memory_limit_mib'}
    base = {k:v for k,v in rows[0]['config'].items() if k not in excluded}
    for study in (prefix, moderate, load(Path('paper-binary/result.json'))):
        for row in study['records']:
            for key in ('vanilla_control', 'cold_start'):
                if key in row:
                    control = row[key]
                    error = max(v['scaled_linf'] for v in control['difference'].values())
                    checks['same_equation_controls'].append(dict(case=study['case'], n=row['n'], control=key,
                        error=error, passed=bool(control['diagnostics']['converged'] and error < 1e-10)))
    checks['physical_free_data_fixed'] = all({k:v for k,v in r['config'].items() if k not in excluded} == base for r in rows)
    hs = [r for r in rows if r['system'] == 'hispid']
    checks['lower_resolution_hispid_trends'] = {region+'_'+norm:bool(all(hs[i+1]['physical'][region][norm] < hs[i]['physical'][region][norm]
        for i in range(len(hs)-1))) for region in ('near', 'bulk') for norm in ('H_rms', 'M_rms')}
    checks['passed_requested_controls'] = bool(checks['physical_free_data_fixed'] and checks['single_hole_analytic_liu']
        and checks['independent_liu_pullback'] and checks['native_single_hole'] and all(x['passed'] for x in checks['same_equation_controls']))
    checks['scope'] = 'solution equivalence, low-grid trends and analytic single-hole validation; no replacement of failed high-grid binary gates'
    checks['input_sha256'] = {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob('*.json') if p.name != 'assessment.json'}
    (root/'assessment.json').write_text(json.dumps(checks, indent=2)+'\n')
    plt.rcParams.update({'font.size':10, 'axes.grid':True, 'grid.alpha':.2})
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    labels = dict(bowen_york='Bowen–York', hispid='Conformal Kerr (QI-wide)')
    for system, color in (('bowen_york', '#275d9c'), ('hispid', '#b84a30')):
        selected = [r for r in rows if r['system'] == system]
        ns = [r['n'][0] for r in selected]
        for ax, region, norm, title in zip(axes.ravel(), ('near', 'bulk', 'near', 'bulk'), ('H_rms', 'H_rms', 'M_rms', 'M_rms'),
                ('Near-hole Hamiltonian', 'Bulk Hamiltonian', 'Near-hole momentum', 'Bulk momentum')):
            ax.loglog(ns, [r['physical'][region][norm] for r in selected], 'o-', color=color, label=labels[system])
            ax.set(title=title, xlabel='Radial count N (polar 2N, azimuthal 16)', ylabel='Physical RMS')
    axes[0, 0].legend(fontsize=9)
    fig.suptitle('Fixed moderate binary: lower grids converge; BY reaches the FD6 observer floor')
    fig.savefig(root/'binary-constraints.png', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), constrained_layout=True)
    for ax, name in zip(axes, ('gamma', 'Kij')):
        for system, color in (('bowen_york', '#275d9c'), ('hispid', '#b84a30')):
            selected = [r for r in rows if r['system'] == system]
            ax.loglog([r['n'][0] for r in selected], [r['difference_from_existing_fine'][name]['scaled_linf'] for r in selected],
                'o-', color=color, label=labels[system])
        ax.set(xlabel='Radial count N', ylabel='max |difference| / (1 + |reference|)', title=name+' versus retained N=224 data')
    axes[0].legend(fontsize=9)
    fig.savefig(root/'binary-solutions.png', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), constrained_layout=True)
    for row in single['single_puncture']:
        axes[0].loglog([x['n'] for x in row['levels']], [x['relative_U_linf'] for x in row['levels']], 'o-',
                      label='BY single '+row['case'])
    for system in ('bowen_york', 'hispid'):
        selected = [r for r in paper['records'] if r['system'] == system]
        axes[1].loglog([r['n'] for r in selected], [r['scalar_relative_U_linf'] for r in selected], 'o-', label=labels[system])
    for ax in axes:
        ax.set(xlabel='Radial count N', ylabel='max relative regularized-U difference')
        ax.legend(fontsize=9)
    axes[0].set_title('Single hole: references N=50 / 70')
    axes[1].set_title('Binary: diagnostic reference N=72')
    fig.savefig(root/'paper-convergence.png', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), constrained_layout=True)
    ns = [r['n'][0] for r in single['liu_extreme']]
    for name in ('gamma', 'Kij'):
        axes[0].semilogy(ns, [r['solution_error'][name] for r in single['liu_extreme']], 'o-', label=name)
    for name in ('H_rms', 'M_rms'):
        line, = axes[1].semilogy(ns, [r['physical'][name] for r in single['liu_extreme']], 'o-', label=name)
        axes[1].axhline(single['liu_exact_floor'][name], linestyle=':', color=line.get_color(), alpha=.7)
    axes[0].set(xlabel='Radial count N (32 polar points)', ylabel='Scaled error versus exact Kerr', title='Liu Kerr χ=0.9999: tensor representation')
    axes[1].set(xlabel='Radial count N (32 polar points)', ylabel='Physical RMS; dotted: exact-seed FD6 floor', title='Independent Cartesian constraints')
    for ax in axes:
        ax.legend(fontsize=9)
    fig.savefig(root/'liu-9999-convergence.png', dpi=180)
    plt.close(fig)
    if a.latex_figures:
        a.latex_figures.mkdir(parents=True, exist_ok=True)
        for name in ('binary-solutions.png', 'liu-9999-convergence.png'):
            shutil.copyfile(root/name, a.latex_figures/name)
    if not checks['passed_requested_controls']:
        raise RuntimeError('a requested equivalence or analytic control failed')


if __name__ == '__main__':
    main()
