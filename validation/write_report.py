"""Generate a compact review report from retained numerical JSON evidence."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def number(value):return f'{value:.6g}'

def main():
    data=json.loads((ROOT/'validation/results.json').read_text())
    lines=['# HiSpID validation report','',
           'All numerical jobs use one CPU thread. Configurations, verifier steps, tensor samples, iterations and runtimes are retained in `validation/results.json` and ignored `validation/raw/`. The native BY baseline is unchanged.',
           '', 'The independent verifier consumes physical gamma/K only. Reported H and the physical momentum norm are in total seed rest-mass units. All binary examples here have total seed rest mass 1. Residual maxima and normalized ratios remain in JSON; normalization denominator floors are 1e-8.',
           '', '## Seed controls','', '| Seed | H RMS | M RMS | Maximum charge error | Passed |', '|---|---:|---:|---:|---|']
    for rec in data.get('seeds',{}).get('records',[]):
        n=rec['step_sequence'][-1]['norms']
        lines.append(f"| {rec['case']} | {number(n['H_rms'])} | {number(n['M_rms'])} | {number(max(rec['charge_error']))} | {rec['passed']} |")
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
        if 'config' in last:
            lines+=['','Fully specified finest-grid input:','', '```json',json.dumps(last['config'],indent=2),'```',
                    '',f"Charge radii: {last.get('charge_radii')}; quadratic inverse-radius extrapolation. Finest [E,P,J]: `{last.get('charges_extrapolated')}`."]
    if 'covariance' in data:
        result=data['covariance'];lines+=['','## Solved coordinate covariance','',f"Gate: **{result.get('passed',False)}**. Source case: `{result.get('source_case','pending')}`."]
        lines+=['','| Grid | gamma relative max | K relative max | Maximum EPJ error |','|---|---:|---:|---:|']
        for rec in result.get('records',[]):
            lines.append('| '+'×'.join(map(str,rec['resolution']))+' | '+number(rec['relative_errors']['gamma'])+' | '+number(rec['relative_errors']['Kij'])+' | '+number(max(rec['charge_error']+rec['origin_charge_error']))+' |')
        lines+=['',f"Translation errors: `{result.get('translation_errors','pending')}`. Rotation errors are compared against five times the measured base-grid truncation difference, with a 1e-9 floor. ADM angular momentum is checked about both translated and fixed origins."]
    lines+=['','## Failures, scope and limitations','',
            '- Initial far40 12²×8,20²×12,28²×16 grids failed badly despite tiny weighted residuals. `validation/failed_moderate_coarse.json` preserves the evidence. The Nyquist derivative defect was independently identified and repaired. Analytic scalar far-filter reparameterization removes the known unresolved shell without changing the equations.',
            '- `validation/failed_seed_verifier_step.json` preserves a verifier-step failure on inner-sheet points; refining the independent stencil restored fourth-order convergence.',
            '- The declared preliminary moderate threshold is RMS 1e-4 with three converging fine grids and <0.5% charge stability. The stronger near/bulk threshold is RMS 1e-6 and max 1e-4. They are reported separately. Neither certifies production or publication accuracy.',
            '- There is no apparent-horizon finder, horizon mass/spin measurement, or binary enclosure certificate. Isolated seed horizon radii are screening information only. Modified regions cannot be called horizon-contained vacuum.',
            '- Exact punctures and signed lapse magnitude below 1e-12 are excluded; extreme-throat robustness and arbitrary-precision arithmetic are not implemented.',
            '- Thesis historical step stuffing differs from modern smooth Eq.26, and boosted thesis spin conventions differ from this rest-spin API. Published high-boost tables omit complete bare inputs. Fully specified local benchmarks must not be called exact table reproduction.',
            '', 'See `docs/HISPID.md` for build/API/conventions, `docs/hispid-review.md` for independent source review, and `docs/HISPID_STATUS.md` for milestone status. Reviewable branch work stays isolated; no main-project merge or shared installation change is performed.', '']
    (ROOT/'docs/VALIDATION.md').write_text('\n'.join(lines))

if __name__=='__main__':main()
