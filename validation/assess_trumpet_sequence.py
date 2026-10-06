"""Assess a declared three-grid sequence from retained independent observers."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
from physical import norms

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def assess(runs,output):
    if len(runs)!=3 or output.exists():raise ValueError('three runs and fresh output required')
    records=[];reference=None;xyz=None
    numerical={'n','tolerance','max_newton','max_krylov','krylov_restart','memory_limit_mib'}
    for path in runs:
        d=json.loads((path/'result.json').read_text());physical={k:v for k,v in d['config'].items() if k not in numerical}
        if reference is None:reference=physical;source=d['library_sha256']
        if reference!=physical or d['library_sha256']!=source:raise ValueError('sequence free data or native image differ')
        with np.load(path/'physical.npz') as r:
            if xyz is None:xyz=r['xyz'].copy()
            if not np.array_equal(xyz,r['xyz']):raise ValueError('physical observer positions differ')
            n=d['physical']['near']['count'];b=d['physical']['bulk']['count'];N=len(xyz)
            bins={'near':np.arange(N)<n,'bulk':(np.arange(N)>=n)&(np.arange(N)<n+b)}
            checked={name:norms(r,mask) for name,mask in bins.items()}
            for name in bins:
                for key in ('H_rms','H_max','M_rms','M_max'):
                    if not np.isclose(checked[name][key],d['physical'][name][key],rtol=1e-12,atol=1e-15):raise ValueError('stored norms disagree with raw observer')
            exterior=bool(np.all(r['stencil_attenuation_all_one'][:n+b]))
        checkpoint=Path(d['checkpoint']['path'])
        if sha(checkpoint)!=d['checkpoint']['file_sha256']:raise ValueError('checkpoint changed')
        passed=bool(d['completed'] and d['diagnostics']['converged'] and d['minimum_psi']>0 and d['minimum_metric_eigenvalue']>0 and exterior and all(v['H_rms']<1e-6 and v['M_rms']<1e-6 and v['H_max']<1e-4 and v['M_max']<1e-4 for v in checked.values()))
        records.append(dict(path=str(path.resolve()),resolution=d['config']['n'],physical=checked,physical_bounds_passed=passed,charges=d['charges_extrapolated'],checkpoint=d['checkpoint'],result_sha256=sha(path/'result.json'),observer_sha256=sha(path/'physical.npz')))
    shapes=np.array([r['resolution'] for r in records]);increasing=bool(np.all(np.diff(shapes,axis=0)>=0) and np.all(np.any(np.diff(shapes,axis=0)>0,axis=1)))
    changes={region:{k:np.diff([r['physical'][region][k] for r in records]).tolist() for k in ('H_rms','M_rms')} for region in ('near','bulk')}
    decreasing=all(all(x<0 for x in v) for row in changes.values() for v in row.values())
    q=np.array([r['charges'] for r in records]);charge_change=np.max(abs(np.diff(q,axis=0))/(1+abs(q[1:])),axis=1)
    out=dict(passed=bool(increasing and decreasing and all(r['physical_bounds_passed'] for r in records) and np.max(charge_change)<1e-3),scope='three-grid physical constraints and charge stability; horizon enclosure separate',binary_acceptance=False,library_sha256=source,records=records,increasing_resolutions=increasing,decreasing_rms=decreasing,rms_changes=changes,charge_scaled_changes=charge_change.tolist(),driver_sha256=sha(Path(__file__)))
    output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));return 0 if out['passed'] else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--runs',type=Path,nargs=3,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();raise SystemExit(assess(a.runs,a.output))
