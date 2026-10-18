"""Analyze a retained Arnoldi projection; no elliptic solve or acceptance."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--plot',type=Path,required=True);a=p.parse_args()
r=json.loads(a.input.read_text().splitlines()[-1]);h=np.asarray(r['hessenberg']);n=h.shape[1];b=np.zeros(n+1);b[0]=1
with np.errstate(all='ignore'):
 q,t=np.linalg.qr(h,mode='reduced');x=np.linalg.solve(t,q.T@b)
 # The generic long-double product avoids a BLAS cancellation diagnostic.
 qr=float(np.linalg.norm(b-h.astype(np.longdouble)@x.astype(np.longdouble)))
 sv=np.linalg.svd(h,compute_uv=False);ritz,V=np.linalg.eig(h[:n]);norm=np.linalg.norm(h[:n])
assert np.isfinite(h).all() and np.isfinite(sv).all() and np.isfinite(qr)
native=r['recurrence_residual']/r['initial_residual'];assert abs(qr-native)<1e-8
out=dict(columns=n,native_recurrence_relative=native,independent_untruncated_qr_relative=qr,
 singular_values=sv.tolist(),projected_operator_norm=float(norm),singular_value_spread=float(sv[0]/sv[-1]),
 eigenvector_condition=float(np.linalg.cond(V)),henrici_departure=float(np.sqrt(max(0,norm*norm-np.sum(abs(ritz)**2)))/norm),
 ritz_real=ritz.real.tolist(),ritz_imag=ritz.imag.tolist(),
 interpretation='Projection is near floating-point conditioning limits and strongly nonnormal. Individual tail singular values/Ritz values are not reliable full-operator estimates. Larger restart reaches the prior residual plateau in one cycle.')
a.output.write_text(json.dumps(out,indent=2)+'\n')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,2,figsize=(10.8,4.3),layout='constrained')
ax[0].semilogy(np.arange(1,n+1),sv,lw=1.8,color='#245c8e')
ax[0].axhline(np.finfo(float).eps*sv[0],color='#ab4539',ls='--',lw=1,label=r'$\epsilon\|\bar H\|_2$')
ax[0].set(xlabel='Singular-value index',ylabel=r'$\sigma_j(\bar H)$',title='Spread approaches floating-point limits');ax[0].legend(frameon=False,fontsize=9);ax[0].grid(alpha=.2)
image=ax[1].imshow(np.ma.masked_less_equal(abs(h),0),origin='upper',aspect='auto',cmap='viridis',norm=LogNorm(vmin=1e-7,vmax=1e13))
ax[1].set(xlabel='Arnoldi column',ylabel='Arnoldi row',title='Large off-diagonal couplings')
fig.colorbar(image,ax=ax[1],label=r'$|\bar H_{ij}|$')
fig.suptitle(r'Original right preconditioner: $240\times480\times16$, one 128-step cycle',fontsize=12)
fig.savefig(a.plot,dpi=180);plt.close(fig)
print({k:v for k,v in out.items() if not isinstance(v,list)})
