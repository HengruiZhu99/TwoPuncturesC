"""Sign audit of the existing interior centered preconditioner stencil."""
import json
from pathlib import Path
import numpy as np
rows=[]
for n,r in ((224,0),(240,0),(224,1),(240,1),(224,4),(240,4)):
 for label,N in (('radial',n),('polar',2*n)):
  h=np.pi/N;angle=h*(np.arange(N)+.5);sn=np.sin(angle);cs=np.cos(angle)
  if label=='radial':
   lam=.2;sigma=(1-cs)/2;d=1-(1-lam)*sigma;t=lam*sigma/d
   dt=lam/d**2;ddt=2*(1-lam)*lam/d**3
   ta=.5*sn*dt;taa=.5*cs*dt+.25*sn**2*ddt
   ctt=t*(1-t)**2;ct=(1-t)*((r+1)*(1-t)-2*t)
   diffusion=ctt/ta**2;drift=ct/ta-ctt*taa/ta**3
  else:
   kap=2.;th=np.tanh(-kap*cs);eta=th/np.tanh(kap)
   de=kap*(1-th**2)/np.tanh(kap);dde=-2*kap**2*th*(1-th**2)/np.tanh(kap)
   eb=de*sn;ebb=de*cs+dde*sn**2
   diffusion=(1-eta**2)/eb**2;drift=-2*(r+1)*eta/eb-(1-eta**2)*ebb/eb**3
  lo=diffusion/h**2-drift/(2*h);hi=diffusion/h**2+drift/(2*h)
  # Both neighbors are distinct, non-reflected nodes on these rows.
  ix=np.flatnonzero((lo[1:-1]<0)|(hi[1:-1]<0))+1
  rows.append(dict(n=N,axis=label,exponent=r,interior_negative_neighbor_rows=ix.tolist(),
     maximum_cell_peclet=float(np.max(abs(drift[1:-1])*h/(2*diffusion[1:-1])))))
Path(__file__).with_suffix('.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
