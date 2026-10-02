"""Generate a compact review report from retained numerical JSON evidence."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def number(value):return f'{value:.6g}'

def main():
    data=json.loads((ROOT/'validation/results.json').read_text())
    lines=['# HiSpID validation report','',
           'All numerical jobs use one CPU thread. Configurations, verifier steps, tensor samples, iterations and runtimes are retained in `validation/results.json` and ignored `validation/raw/`. The native BY baseline is unchanged.',
           '', 'The independent verifier differentiates physical gamma/K only; g is metadata for the region bins. Reported H and the physical momentum norm are in total seed rest-mass units. All binary examples here have total seed rest mass 1. RMS means an average over fixed off-grid points, not a volume L2 norm. Residual maxima and normalized ratios remain in JSON; normalization denominator floors are 1e-8. The normalized momentum ratio is uninformative for maximal data, so absolute physical norms define acceptance.',
           '', 'Numerical evidence is tied to each record\'s native library SHA. Older failed-source cases remain visible. API migration compares fresh, separate native processes with loaded-image checks; checkpoint guards otherwise remain strict. The original same-process comparison was invalid because dyld reused an archived image with an identical install name. It and the affected attempted revalidation are preserved and explicitly withdrawn.',
           '', 'Current Cartesian-regular modal-P binaries remain unaccepted. The leading current source uses mapped Chebyshev coordinates and exact modal FD block elimination; old nodal-V binary results are historical and do not establish axis regularity. All tables retain their own source fingerprints and failed flags.',
           '', '## Seed controls','', '| Seed | H RMS | M RMS | Maximum charge error | Passed |', '|---|---:|---:|---:|---|']
    for rec in data.get('seeds',{}).get('records',[]):
        n=rec['step_sequence'][-1]['norms']
        lines.append(f"| {rec['case']} | {number(n['H_rms'])} | {number(n['M_rms'])} | {number(max(rec['charge_error']))} | {rec['passed']} |")
    target=ROOT/'validation/target_seed_controls_current.json'
    if not target.exists():target=ROOT/'validation/target_seed_controls.json'
    if target.exists():
        lines+=['','## Revised isolated seed targets','',
                'Seed rest spin chi=.95 and lab speed v=.885, separately and together with generic directions. This does not accept a solved binary. Independent Cartesian constraints use three verifier step sizes; charge fits use increasing extraction radii with separately refined angular quadrature.',
                '', '| Case | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | Maximum charge error | Passed |',
                '|---|---:|---:|---:|---:|---:|---|']
        for rec in json.loads(target.read_text())['cases']:
            last=rec['constraint_sequence'][-1]
            values=[last[k][q] for k in ('near','bulk') for q in ('H_rms','M_rms')]+[rec['charge_error_max']]
            lines.append('| '+rec['case']+' | '+' | '.join(map(number,values))+f" | {rec['passed']} |")
    for label,result in data.items():
        if label in ('seeds','metadata','covariance') or not isinstance(result,dict) or not result.get('records'):continue
        lines+=['',f'## {label}','',f"Preliminary gate: **{result.get('passed',False)}**. Stronger gate: **{result.get('passed_strict',False)}**. Horizon enclosure: **unverified**.",
                '', '| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |', '|---|---:|---:|---:|---:|---:|---:|']
        for rec in result['records']:
            if 'near' not in rec:continue
            grid='×'.join(map(str,rec['resolution']))
            values=[rec['near']['H_rms'],rec['near']['M_rms'],rec['bulk']['H_rms'],rec['bulk']['M_rms'],rec['charges_extrapolated'][0],rec['diagnostics']['seconds']]
            lines.append('| '+grid+' | '+' | '.join(map(number,values))+' |')
        lines+=['','| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |','|---|---:|---:|---:|---:|---|']
        for rec in result['records']:
            if 'near' not in rec:continue
            grid='×'.join(map(str,rec['resolution']));atten=rec['attenuation'];ext=rec['g_equals_one']
            inside=[number(atten[k]) if atten.get('count',0) else 'no support' for k in ('H_rms','M_rms')]
            internal=', '.join(map(number,rec.get('physical_equivalent_g1_linf',[]))) or 'not captured'
            lines.append('| '+grid+' | '+' | '.join(inside+[number(ext['H_max']),number(ext['M_max']),internal])+' |')
        last=result['records'][-1]
        if 'library_sha256' in last:
            lines+=['',f"Finest source SHA: `{last['library_sha256']}`. Unknown basis: `{last.get('unknown_parameterization_id',last.get('unknown_parameterization','historical'))}`; maps: `{last.get('collocation_maps')}`."]
        if result.get('reference_comparison'):
            lines+=['', 'Reference comparison:', '', '```json',json.dumps(result['reference_comparison'],indent=2),'```']
        if 'accepted_high_regime' in result:
            lines+=['',f"Combined high-regime gate: **{result['accepted_high_regime']}**. Finest local strict gate: **{last.get('passed_strict',False)}**. Original three-grid monotonic gate: **{result.get('converges',False)}**."]
        if 'config' in last:
            lines+=['','Fully specified finest-grid input:','', '```json',json.dumps(last['config'],indent=2),'```',
                    '',f"Charge radii: {last.get('charge_radii')}; quadratic inverse-radius extrapolation. Finest [E,P,J]: `{last.get('charges_extrapolated')}`."]
    if 'covariance' in data:
        result=data['covariance'];lines+=['','## Solved coordinate covariance','',f"Gate: **{result.get('passed',False)}**. Source case: `{result.get('source_case','pending')}`."]
        lines+=['','| Grid | gamma relative max | K relative max | Maximum EPJ error |','|---|---:|---:|---:|']
        for rec in result.get('records',[]):
            lines.append('| '+'×'.join(map(str,rec['resolution']))+' | '+number(rec['relative_errors']['gamma'])+' | '+number(rec['relative_errors']['Kij'])+' | '+number(max(rec['charge_error']+rec['origin_charge_error']))+' |')
        lines+=['',f"Translation errors: `{result.get('translation_errors','pending')}`. Rotation errors are compared against five times the measured base-grid truncation difference, with a 1e-9 floor. ADM angular momentum is checked about both translated and fixed origins. The 12×24 surface rule has a rotation-dependent angular quadrature error; see the separate refined-charge evidence below."]
    floor_path=ROOT/'validation/verifier_floor_highspin.json'
    if floor_path.exists():
        floor=json.loads(floor_path.read_text())
        lines+=['','## High-spin verifier calibration','',
                'The original all-norms monotonic gate remains failed. On the identical 18 bulk points, exact vacuum Brill–Lindquist and both chi=.99 Kerr seeds show the same worsening with smaller Cartesian stencils as the solved binary. This supports a numerical-floor interpretation of the bulk plateau; it does not measure the true bulk residual or change the acceptance gate.',
                '', '| Exact control or solved grid | h=.016 H RMS | h=.008 | h=.004 | h=.002 | h=.001 |','|---|---:|---:|---:|---:|---:|']
        for rec in floor['controls']+floor['binary_grids']:
            label=rec.get('case',str(rec.get('resolution')))
            lines.append('| '+label+' | '+' | '.join(number(r['norms']['H_rms']) for r in rec['sequence'])+' |')
        lines+=['','Full pointwise/stencil evidence: `validation/verifier_floor_highspin.json`. High-spin angular refinement from160²×16 to160²×24 lowers near H RMS from1.55e-6 to3.18e-7 and M RMS from8.69e-8 to2.70e-9. Near stencil checks leave those results stable. Local strict thresholds pass, but the predeclared aggregate gate does not.']
    lines+=['','## Current integration and conditioning diagnostics','',
            'The matched default-map80×80→80×160×28 polar refinement reduces fixed bulk physical momentum RMS5.79x and identical104×208×28 g1 momentum RMS5.51x. Both binaries remain failed. The new80×160/104×208/128×256×28 sequence is declared in polar_sequence_plan.json; its current records are shown above. Physical thresholds are unchanged.',
            '', 'Analytic metric gradients reduce matched charge integration from13.657s to1.120s (12.20x), with energy differences<2.7e-11 and bit-identical P/J on matching sphere nodes. Global-z polar refinement at fixed phi64 gave misleading stability: phi128 shifted angular momentum by.00935. Aligning the sphere polar axis with the prolate axis yields azimuthal change6.19e-14 and polar change6.35e-10 on the same failed polynomial. Matched independent Cartesian-FD fluxes and centered tensor covariance pass1e-8, with discrepancies<5.6e-11 and1.7e-10. These are integration controls, not binary acceptance.',
            '', 'The optional sin6/(1-t)^6 positive row scaling fails its declared exact-seed far-source floor1e-14: the finest boosted/generic controls reach1.89e-13. No binary has been solved in that norm. The failed artifact is retained. Correcting the stale76-point allocation estimate to the current five-point stencil preserves fields/residual/JVP bit-for-bit and keeps the2048MiB default and8192MiB cap.',
            '', 'The three-grid moderate polar sequence subsequently passes the preliminary physical, refined centered charge and solved coordinate-covariance gates. Current moderate direct AthenaK horizons and inner-ball enclosure also pass; stronger accuracy and revised high-parameter binaries remain unverified. Centered extraction uses2Npolar/3Npolar and phi64/128; angular/fit changes are<1.8e-13, and the two finest extrapolated EPJ vectors differ2.33e-7. Ring reuse gives another21.12x speedup on matched64×64 extraction, with differences<5.6e-13; all-mode first-gradient, frame, invalidation and independent centered/off-center FD controls pass.',
            '', 'The covariance fixed-global-origin charges still use coarse12×24 integration: finest error1.99e-4 passes the preliminary.005 bound and is not a stronger quadrature claim. Its process peak is10.32GB over sequential contexts. See docs/HISPID_PERFORMANCE.md for measured serial speedups and the isolated integration handoff.',
            '', '## Additional reproducible checks','']
    for name in ('asymptotic_inverse_default.json','analytic_charge_method_comparison.json','analytic_charges_operator_equivalence.json','aligned_charges_operator_equivalence.json','aligned_charge_fd_comparison.json','budget_operator_equivalence.json','charge_quadrature_polar2_r128.json','charge_quadrature_polar2_r128_fine.json','charge_quadrature_polar2_r128_analytic.json','charge_quadrature_polar2_r128_phi_refinement.json','charge_quadrature_polar2_aligned.json','far_source_floor_infinity_equilibrated.json','polar_sequence_plan.json','charge_rings_operator_equivalence.json','charge_rings_method_comparison.json','charge_rings_independent_fd.json','charge_rings_private_controls.json','charge_rings_full_tests.json','polar_sequence_refined_charges.json','polar_sequence_horizons.json','preconditioner_reuse_initial_fixture_failure.json','preconditioner_reuse_full_tests.json','preconditioner_reuse_operator_equivalence.json','preconditioner_reuse_benchmark.json'):
        path=ROOT/'validation'/name
        if path.exists():
            value=json.loads(path.read_text());flag=value.get('passed',value.get('passed_off_axis_and_equations','diagnostic only'))
            lines.append(f'- `validation/{name}`: control flag={flag}; binary acceptance is not implied.')
    for name in ('highboost_seed_controls.json','failed_highboost_seed_charge_resolution.json','target_seed_controls.json','target_seed_controls_current.json','target_seed_controls_failed_short_radial_fit.json','charge_quadrature_moderate.json','charge_quadrature_highspin.json','highboost_continued_far0_continuation.json','highboost_continued_fine_far0_continuation.json','highboost_actualop_far0_continuation.json','axis_api_migration.json','failed_axis_basis_trials.json','regular_modes_mapped.json','regular_operators_difference.json','mapped_cache_equivalence.json','modal_block_equivalence.json','difference_derivative_equivalence.json','focused_map_operator_audit.json','private_operator_controls_focus05_k3.json','regular_operators_focus05_k3.json','equation_rows_focus05_k3.json','dense_equations_focus05_k3.json','dense_equations_phi_focus05_k3.json','dense_equations_angular_focus05_k3.json','dense_equations_radial_focus05_k3.json','collocation_map_api_equivalence.json','anisotropic_inverse_default.json','anisotropic_inverse_restart128.json','dense_equations_default_square_r128.json','dense_equations_default_polar2_r128.json','angular_refinement_comparison.json'):
        path=ROOT/'validation'/name
        if path.exists():
            value=json.loads(path.read_text())
            if 'stages' in value:
                description='; '.join(f"v={r['velocity']:.6g}: status={r['diagnostics']['status']}, Newton={r['diagnostics']['newton_iterations']}, GMRES={r['diagnostics']['krylov_iterations']}" for r in value['stages'])
            elif 'passed_off_axis_and_equations' in value:
                description=f"off-axis/equation compatibility={value['passed_off_axis_and_equations']}; solved-axis stencil limitations are recorded separately"
            elif 'passed' in value:description=f"passed={value['passed']}"
            else:description='independent and refined angular/radial charge evidence'
            lines.append(f'- `validation/{name}`: {description}.')
    lines+=['','## AthenaK integration','',
            'The isolated AthenaK branch starts at PR790 head22baa243970fa1880b2bbc48e88a590069d55e47. Its z4c/hispid pgen loads portable physical gamma/K checkpoints and checks ADM/Z4c round trips. AthenaK FastFlow finds exact isolated Schwarzschild, chi=.95 Kerr, v=.885 boosted Schwarzschild and combined chi=.95/v=.885 Kerr horizons at initial time with zero evolution steps.',
            '', 'The combined seed uses flow alpha=.2,lmax48,ntheta50: expansion RMS9.92e-8, relative area error1.75e-14 and independent sampled relative shape error3.35e-7. Fixed-order ntheta74 quadrature confirms area/RMS stability. The failed alpha1 run and coarse strict flags are retained. The three-level outer mesh constraint checks pass at64³, with combined H/M RMS6.13e-7/2.23e-7. These tests validate isolated dataset import and finder controls, not binary attenuation enclosure.',
            '', 'Current executable/source fingerprints and replay records are in the sibling AthenaK worktree docs/hispid-current-controls.json (native9cbf1108…); historical sampler controls and fixed-order quadrature are in docs/hispid-validation.json. Build and replay instructions are in its docs/hispid.md.']
    lines+=['','## Failures, scope and limitations','',
            '- Initial far40 12²×8,20²×12,28²×16 grids failed badly despite tiny weighted residuals. `validation/failed_moderate_coarse.json` preserves the evidence. The Nyquist derivative defect was independently identified and repaired. Analytic scalar far-filter reparameterization removes the known unresolved shell without changing the equations.',
            '- `validation/failed_seed_verifier_step.json` preserves a verifier-step failure on inner-sheet points; refining the independent stencil restored fourth-order convergence.',
            '- The far40 moderate source remains unsupported by the three-grid physical gate (`failed_moderate_filtered_after_split.json`). The original high-spin source suffered singular seed-divergence cancellation; failed H/M trends are retained. The analytic seed momentum identity repaired that source and restored symmetric charges and decreasing momentum norms.',
            '- Direct Gamma=sqrt5 binary solves failed the nonlinear/Krylov and physical gates; `highboost_far0` retains all three resolutions. Continuation records distinguish solver convergence from independent physical validation. No boost interval is inferred from a few successful Newton stages.',
            '- The declared preliminary moderate threshold is RMS 1e-4 with three converging fine grids and <0.5% charge stability. The stronger near/bulk threshold is RMS 1e-6 and max 1e-4. They are reported separately. Neither certifies production or publication accuracy.',
            '- Exact isolated seed and current moderate-binary direct horizon controls pass. Moderate inner g/operator balls are enclosed by continuous retained-surface bounds with empirical refinement allowance. Generic AKV spin, stronger exterior accuracy and high-parameter binary horizons remain unverified.',
            '- Historical unconstrained nodal-V fields violate Cartesian axis regularity. Current mapped modal-P fields enforce C2 axis limits and pass independent scalar/vector Cartesian manufactured checks. The current polar sequence passes preliminary physical/charge/covariance/horizon gates. Its stronger physical gate still fails; regularity and small internal residuals do not establish production exterior accuracy.',
            '- Exact punctures remain excluded. Regular signed-lapse graph formulas support the QI throat and pass independent seed/derivative checks there; high-regime binary accuracy and arbitrary-precision arithmetic remain separate questions.',
            '- Thesis historical step stuffing differs from modern smooth Eq.26, and boosted thesis spin conventions differ from this rest-spin API. Published high-boost head-on descriptions omit bare masses and companion-attenuation widths. No original parameter files were recovered. Fully specified local benchmarks must not be called exact table reproduction; milestone E remains incomplete.',
            '', 'See `docs/HISPID.md` for build/API/conventions, `docs/hispid-review.md` for independent source review, and `docs/HISPID_STATUS.md` for milestone status. Reviewable branch work stays isolated; no main-project merge or shared installation change is performed.', '']
    (ROOT/'docs/VALIDATION.md').write_text('\n'.join(lines))

if __name__=='__main__':main()
