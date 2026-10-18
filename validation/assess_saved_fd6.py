"""Assess the declared FD6 observers without altering the retained FD4 verdict."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from physical import norms
p=argparse.ArgumentParser();p.add_argument('--runs',type=Path,required=True);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
records=[];fields=[];original=[];same=None;image=None;pointset=None
for n in (160,192,224):
 root=a.runs/f'hispid{n}';r=json.loads((root/'result.json').read_text());s=json.loads((a.source/f'moderate{n}/result.json').read_text());cfg=s['config'].copy();cfg.pop('n')
 if same is None:same=cfg;image=r['library_sha256']
 assert cfg==same and r['library_sha256']==image and s['completed'] and s['diagnostics']['converged']
 with np.load(root/'physical.npz') as q:
  if pointset is None:pointset=q['xyz'].copy()
  assert np.array_equal(pointset,q['xyz']) and np.all(q['stencil_attenuation_all_one']) and np.min(q['min_metric_eigenvalue'])>0
  checked={name:norms(q,mask) for name,mask in [('near',np.arange(66)<48),('bulk',np.arange(66)>=48)]}
  for name in checked:
   for key in ('H_rms','M_rms','H_max','M_max'):assert np.isclose(checked[name][key],r['physical'][name][key],rtol=1e-12,atol=1e-15)
  fields.append((q['gamma'].copy(),q['Kij'].copy()))
 records.append(dict(resolution=r['resolution'],physical=checked,physical_sha256=r['physical_sha256'],source_result_sha256=r['source_result_sha256']))
 original.append(s)
changes={region:{key:np.diff([r['physical'][region][key] for r in records]).tolist() for key in ('H_rms','M_rms')} for region in ('near','bulk')}
decrease=all(all(d<0 for d in values) for row in changes.values() for values in row.values())
bounds=all(q['H_rms']<1e-6 and q['M_rms']<1e-6 and q['H_max']<1e-4 and q['M_max']<1e-4 for row in records for q in row['physical'].values())
charges=np.array([s['charges_extrapolated'] for s in original]);delta=np.max(abs(np.diff(charges,axis=0))/(1+abs(charges[1:])),axis=1)
field_changes={name:[float(np.max(abs(fine[k]-coarse[k])/(1+abs(fine[k])))) for coarse,fine in zip(fields,fields[1:])] for k,name in enumerate(('gamma','Kij'))}
fd4=json.loads((a.source/'assessment.json').read_text());calibration=json.loads((a.runs/'control.json').read_text());assert calibration['passed']
result=dict(passed=bool(bounds and decrease and np.max(delta)<1e-3),scope='Fixed physical formulation and three saved grids; predeclared higher-accuracy Cartesian observer. Original FD4 assessment retained.',binary_acceptance=False,fd4_passed=fd4['passed'],physical_bounds_passed=bounds,decreasing_rms=decrease,rms_changes=changes,charge_scaled_changes=delta.tolist(),field_cauchy_changes=field_changes,records=records,library_sha256=image,calibration_sha256=hashlib.sha256((a.runs/'control.json').read_bytes()).hexdigest())
a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));raise SystemExit(0 if result['passed'] else 1)
