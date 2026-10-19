"""Fourier localization of a retained nonlinear residual; no solve or acceptance transfer."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python')]
from hispid import Backend
from checkpoint_export import read_checkpoint

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--library',type=Path,required=True);p.add_argument('--checkpoint',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();b=Backend(str(a.library.resolve()));c,u,meta=read_checkpoint(a.checkpoint)
if meta['parameterization']!=b.parameterization():raise ValueError('checkpoint/basis mismatch')
with b.create(c,execution='kokkos') as s:
    s.set_unknowns(u);r=s.residual(u).reshape(c.n[2],c.n[1],c.n[0],4)
f=np.fft.rfft(r,axis=0,norm='ortho');weight=np.full(len(f),2.);weight[[0,-1]]=1.
power=weight[:,None]*np.sum(abs(f)**2,axis=(1,2));total=float(np.sum(r*r));rows=[]
for m in range(len(f)):
    j,i,v=np.unravel_index(np.argmax(abs(f[m])),f[m].shape)
    rows.append(dict(m=m,l2_by_component=np.sqrt(power[m]).tolist(),total_l2_fraction=float(np.sum(power[m])/total),maximum_nodal_index=dict(radial=int(i),polar=int(j),component=int(v))))
result=dict(purpose='scaled native residual localization only, not physical constraints or convergence',checkpoint=meta,library_sha256=b.library_sha256(),driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            component_order=['Hamiltonian','local_vector0','local_vector1','local_vector2'],scaled_linf_by_component=np.max(abs(r),axis=(0,1,2)).tolist(),l2=float(np.sqrt(total)),parseval_relative_error=float(abs(np.sum(power)-total)/total),rows=rows)
a.output.write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)
