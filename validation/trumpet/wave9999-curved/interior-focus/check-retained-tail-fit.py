"""Read-only held-out test of a bounded far-zone restart ansatz; no checkpoint edit."""
import argparse,json
from pathlib import Path
import numpy as np
from checkpoint_export import read_checkpoint
a=argparse.ArgumentParser();a.add_argument('checkpoint');a.add_argument('output');args=a.parse_args()
c,p,meta=read_checkpoint(args.checkpoint);na,nb,np_=c.n
assert meta['file_sha256']=='c3132b0510ac766c5ccc31ac0e2a75cedf24ff6103cc9508457cbe4212acda7a', 'This evidence replay requires the exact retained checkpoint'
p=p.reshape(np_,nb,na,4)
sigma=.5*(1-np.cos(np.pi*(np.arange(na)+.5)/na));r=16*sigma/(1-sigma)
fit=np.where((r>=256)&(r<=1024))[0];train=fit[::2];held=fit[1::2]
x=1024/(r+.5);V=np.polynomial.polynomial.polyvander(x,2);rows=[]
for m in range(1,5):
 slots=[m,np_//2+m];data=p[slots,:,:,1:];y=data.transpose(2,0,1,3).reshape(na,-1)
 coeff=np.linalg.lstsq(V[train],y[train],rcond=None)[0];pred=np.einsum("ij,jk->ik",V,coeff,optimize=False)
 assert np.isfinite(pred).all()
 err=np.linalg.norm(pred[held]-y[held])/np.linalg.norm(y[held]);en=np.linalg.norm(pred[held]-y[held],axis=0)
 rows.append(dict(mode=m,held_out_relative=float(err),held_out_maximum_error=float(np.max(abs(pred[held]-y[held]))),held_out_maximum_value=float(np.max(abs(y[held]))),fitted_infinity_maximum=float(np.max(abs(coeff[0]))),actual_outermost_maximum=float(np.max(abs(y[-1]))),passed=bool(err<.01)))
out=dict(scope='Initial-guess tail model diagnostic only; no original checkpoint modification and no solve or acceptance',checkpoint_sha256=meta['file_sha256'],fit_band=[256,1024],degree=2,train_indices=train.tolist(),held_out_indices=held.tolist(),design_condition=float(np.linalg.cond(V[train])),criterion='Every mode held-out relative P error <.01; this is not a physical convergence tolerance',rows=rows,passed=all(v['passed'] for v in rows))
Path(args.output).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
