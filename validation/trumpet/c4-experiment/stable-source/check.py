"""Check source rearrangement against retained equations and a directional derivative."""
import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'examples')]
from hispid import Backend
from trumpet_configs import trumpet_moderate
p=Path(__file__).resolve().parent
old=p.parent/'scaled-inverse'
a=json.loads((old/'baseline-fields.json').read_text());new=json.loads((p/'control.json').read_text())
with np.load(old/'baseline-fields.npz') as x,np.load(p/'control.npz') as y:
 result=dict(residual_linf=float(np.max(abs(x['residual']-y['residual']))),jvp_linf=float(np.max(abs(x['jvp']-y['jvp']))))
result['converged_field_scaled_linf']=max(float(np.max(abs(np.asarray(a['field_witness'][k])-np.asarray(new['field_witness'][k]))/(1+abs(np.asarray(a['field_witness'][k]))))) for k in a['field_witness'])
b=Backend(str(ROOT/'build-c4-stable-source/libHiSpID.dylib'))
c=trumpet_moderate(b,12,16)
with b.create(c) as s:
 k=np.arange(s.size);u=1e-3*np.sin(.13*k);d=1e-3*np.cos(.27*k);j=s.jvp(u,d)
 result['directional_errors']=[]
 for h in (.02,.01,.005):
  fd=(s.residual(u+h*d)-s.residual(u-h*d))/(2*h)
  result['directional_errors'].append(dict(h=h,relative_l2=float(np.linalg.norm(fd-j)/np.linalg.norm(j))))
result['binary_acceptance']=False
assert result['residual_linf']<1e-13 and result['jvp_linf']==0
assert result['directional_errors'][-1]['relative_l2']<1e-8
assert new['diagnostics']['converged']
# Field drift is reported, not hidden behind an equation-level tolerance.
(p/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
