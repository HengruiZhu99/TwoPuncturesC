"""Compare one predeclared polar refinement on identical PDE/Cartesian points.

This postprocessor creates no native context. Auxiliary P angular tails are
reported descriptively; they do not replace physical constraint acceptance.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from run_validation import SOLVER_CONTROLS

ROOT=Path(__file__).resolve().parents[1]

def auxiliary_tails(record):
    na,nb,np_=record['resolution']
    path=ROOT/f"validation/raw/{record['case']}_{na}_{np_}.npz"
    with np.load(path) as saved:values=saved['unknowns'].reshape(np_,nb,na,4)
    degree=np.arange(nb)
    matrix=2/nb*(-1.)**degree[:,None]*np.cos(np.pi*degree[:,None]*(np.arange(nb)+.5)/nb)
    coefficients=np.einsum('ji,kilv->kjlv',matrix,values)
    coefficients[:,0]*=.5
    all_norm=np.sqrt(np.sum(coefficients**2,axis=(0,1,2)))
    tail_norm=np.sqrt(np.sum(coefficients[:,int(np.ceil(.8*nb)):]**2,axis=(0,1,2)))
    return dict(top_degree_fraction=.2,relative_l2_by_variable=(tail_norm/np.maximum(all_norm,1e-300)).tolist(),
                raw_checkpoint_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                note='Raw auxiliary P Chebyshev coefficients; not a physical volume norm or acceptance gate.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--baseline-dense',required=True)
    p.add_argument('--refined-dense',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();plan=json.loads((ROOT/'validation/angular_refinement_plan.json').read_text())
    all_results=json.loads((ROOT/'validation/results.json').read_text())
    records=[all_results[plan[key]]['records'][-1] for key in ('baseline_case','refined_case')]
    dense=[json.loads(Path(path).read_text()) for path in (a.baseline_dense,a.refined_dense)]
    for record,replay in zip(records,dense):
        if record['library_sha256']!=plan['library_sha256'] or replay['library_sha256']!=plan['library_sha256']:
            raise ValueError('native source differs from the predeclared plan')
        if replay['case']!=record['case'] or replay['source_resolution']!=record['resolution']:
            raise ValueError('dense replay differs from selected checkpoint')
        if replay['evaluation_resolution']!=plan['comparison_grid'] or replay['solve_performed']:
            raise ValueError('comparison requires the identical declared grid, without another solve')
    physical_inputs=[{k:v for k,v in record['config'].items() if k not in SOLVER_CONTROLS} for record in records]
    if physical_inputs[0]!=physical_inputs[1]:raise ValueError('free data or equation choices differ')
    if dense[0]['groups']['g1']['count']!=dense[1]['groups']['g1']['count']:
        raise ValueError('exterior sample masks differ')
    native_momentum=[float(np.linalg.norm(replay['groups']['g1']['physical_equivalent_rms'][1:])) for replay in dense]
    bulk_momentum=[record['bulk']['M_rms'] for record in records]
    reductions=[native_momentum[0]/native_momentum[1],bulk_momentum[0]/bulk_momentum[1]]
    converged=all(record['diagnostics']['status']==0 for record in records)
    tails=[auxiliary_tails(record) for record in records]
    growing_tail=any(y>x for x,y in zip(tails[0]['relative_l2_by_variable'],tails[1]['relative_l2_by_variable']))
    conclusion=('failure' if min(reductions)<2 or growing_tail else 'success' if min(reductions)>=4 else 'inconclusive') if converged else 'unconverged'
    output=dict(plan=plan,records=[dict(case=r['case'],resolution=r['resolution'],diagnostics=r['diagnostics'],
        near=r['near'],bulk=r['bulk'],angular_auxiliary_tails=t) for r,t in zip(records,tails)],
        common_grid_g1_cartesian_momentum_components_rms_norm=native_momentum,
        fixed_bulk_physical_momentum_norm_rms=bulk_momentum,
        reduction_factors=dict(common_grid=native_momentum[0]/native_momentum[1],bulk=reductions[1]),
        growing_auxiliary_angular_tail=growing_tail,diagnostic_conclusion=conclusion,acceptance=False,
        note='One directional diagnostic, not a convergence sequence. Physical/charge/axis/covariance/horizon acceptance thresholds are unchanged.')
    Path(a.output).write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))

if __name__=='__main__':main()
