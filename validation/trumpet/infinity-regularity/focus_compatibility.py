"""Read-only low-mode focus regularity diagnostic for the original radial map.

For m<=4, u_m/rho^m = -2(1-t)^(m+1) P_m/(2b)^m. Cartesian
analyticity at a focus implies (d_t-2*side*d_eta-(m+1)) P_m=0.
This evaluates that necessary condition, not a physical constraint norm.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
if a.output.exists():raise FileExistsError(a.output)
r=json.loads((a.run/'result.json').read_text());na,nb,np_=r['config']['n']
with (a.run/'diagnostic.checkpoint').open() as f:
 f.readline();parameterization=f.readline().split()[1]
assert parameterization=='modal_P_C2tauC4_map_v5_r0.20000000000000001_k2'
assert hashlib.sha256((a.run/'solve.npz').read_bytes()).hexdigest()==r['solve_artifact_sha256']
with np.load(a.run/'solve.npz') as z:u=z['unknowns'].reshape(np_,nb,na,4).astype(np.longdouble)
def endpoints(n,side):
 theta=np.arccos(np.longdouble(-1))*(np.arange(n,dtype=np.longdouble)+.5)/n
 delta=-2*np.sin(theta/2)**2 if side<0 else 2*np.cos(theta/2)**2
 w=(-1.)**np.arange(n)*np.sin(theta)/delta;w/=np.sum(w)
 d=w*(np.sum(w/delta)-1/delta)
 return w,d
wa,da=endpoints(na,-1);da*=2/np.longdouble(.2)
rows=[]
for side in (-1,1):
 w,d=endpoints(nb,side);d*=np.tanh(np.longdouble(2))/(2*(1-np.tanh(np.longdouble(2))**2))
 for mode in range(np_):
  m=mode if mode<=np_//2 else mode-np_//2
  if m>4:continue
  value=np.einsum('j,i,jiv->v',w,wa,u[mode]);dt=np.einsum('j,i,jiv->v',w,da,u[mode]);de=np.einsum('j,i,jiv->v',d,wa,u[mode])
  defect=dt-2*side*de-(m+1)*value
  rows.append(dict(side=side,mode=mode,m=m,defect=defect.astype(float).tolist(),relative=(abs(defect)/(1+abs(dt)+2*abs(de)+(m+1)*abs(value))).astype(float).tolist()))
a.output.write_text(json.dumps(dict(scope=__doc__,parameterization=parameterization,rows=rows),indent=2)+'\n')
