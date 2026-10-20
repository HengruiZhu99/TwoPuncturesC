"""Inspect the first, unaugmented vector-inner cycle; no solve or acceptance."""
import argparse,gzip,hashlib,json
from pathlib import Path
import numpy as np
from scipy import linalg
p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
with (gzip.open(a.input,'rb') if a.input.suffix=='.gz' else a.input.open('rb')) as stream:line=next(stream).rstrip(b'\r\n')
r=json.loads(line)
assert r['method']==2 and r['cycle_begin']==0 and 1<=r['columns']<=256
h=np.asarray(r['hessenberg']);n=r['columns'];assert h.shape==(n+1,n) and np.isfinite(h).all()
b=np.zeros(n+1);b[0]=1
q,t=linalg.qr(h,mode='economic');x=linalg.solve_triangular(t,q[0,:])
qr=float(np.linalg.norm(b.astype(np.longdouble)-h.astype(np.longdouble)@x.astype(np.longdouble)))
u,s,vh=linalg.svd(h,full_matrices=True,lapack_driver='gesvd')
columns=linalg.norm(h,axis=0);balanced=h/columns
sb=linalg.svdvals(balanced)
q2,t2=linalg.qr(balanced,mode='economic');y=linalg.solve_triangular(t2,q2[0,:])
balanced_qr=float(np.linalg.norm(b.astype(np.longdouble)-balanced.astype(np.longdouble)@y.astype(np.longdouble)))
eps=np.finfo(float).eps
out=dict(record='first unaugmented vector-inner cycle',raw_record_sha256=hashlib.sha256(line).hexdigest(),columns=n,native_recurrence_relative=r['recurrence_residual']/r['initial_residual'],independent_untruncated_qr_relative=qr,svd_left_null_relative_unreliable=float(abs(u[0,-1])),operator_2norm=float(s[0]),singular_value_spread=float(s[0]/s[-1]),column_norm_min=float(min(columns)),column_norm_max=float(max(columns)),column_norm_spread=float(max(columns)/min(columns)),column_balanced_spread=float(sb[0]/sb[-1]),column_balanced_qr_relative=balanced_qr,singular_values=s.tolist(),normalized_singular_values=(s/s[0]).tolist(),singular_values_below_eps_norm=int(np.sum(s<=eps*s[0])),singular_values_below_n_eps_norm=int(np.sum(s<=n*eps*s[0])),rhs_left_singular_projection_squared=(u[0,:]**2).tolist(),weak_direction_gap_4=float(s[-5]/s[-4]) if n>4 else None,weak4_rhs_projection_squared=float(np.sum(u[0,n-4:n]**2)),interpretation='Small projected-system diagnostic only. Tail singular values near rounding bounds do not establish reliable full-Jacobian eigenmodes. Column balancing leaves the exact search space unchanged. Original native and physical acceptance gates are unchanged.')
a.output.write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if not isinstance(v,list)})
