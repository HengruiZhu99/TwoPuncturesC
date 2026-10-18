"""Read-only C1/C2 Cartesian compatibility of retained modal unknowns."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--runs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();records=[]
R=np.longdouble;pi=np.arccos(R(-1));lam=R('.2');kap=R(2);T=np.tanh(kap)
def weights(n,side):
 th=pi*(np.arange(n,dtype=R)+R('.5'))/n;dx=2*np.sin(th/2)**2 if side<0 else -2*np.cos(th/2)**2
 w=(-1.)**np.arange(n)*np.sin(th)/(-dx);w/=np.sum(w);q=-1/dx;S=np.sum(w*q);S2=np.sum(w*q*q)
 return w,w*(S-q),2*w*(S*S-S*q+q*q-S2)
for n in (160,192,224):
 root=a.runs/f'moderate{n}';r=json.loads((root/'result.json').read_text());assert r['config']['far_radius']==0 and r['config']['inner_flatten']
 assert hashlib.sha256((root/'solve.npz').read_bytes()).hexdigest()==r['solve_artifact_sha256']
 na,nb,np_=r['config']['n'];u=np.load(root/'solve.npz')['unknowns'].reshape(np_,nb,na,4).astype(R)
 wa,az,azz=weights(na,-1);da=2/lam*az;daa=4/(lam*lam)*azz-4*(1-lam)/(lam*lam)*az
 row=[]
 for side in (-1,1):
  wb,bz,bzz=weights(nb,side);s=T/(kap*(1-T*T));db=s*bz;dbb=s*s*bzz+2*side*T**3/(kap*(1-T*T)**2)*bz
  for k in range(np_):
   m=k if k<=np_//2 else k-np_//2
   if m>4:continue
   c1=wb[:,None]*da[None,:]-2*side*db[:,None]*wa[None,:]-(m+1)*wb[:,None]*wa[None,:]
   c2=wb[:,None]*daa[None,:]-4*dbb[:,None]*wa[None,:]-2*(m+1)*wb[:,None]*da[None,:]-4*side*db[:,None]*wa[None,:]+m*(m+1)*wb[:,None]*wa[None,:]
   row.append(dict(side=side,mode=k,m=m,C1=(np.sum(c1[:,:,None]*u[k],axis=(0,1))/np.sum(abs(c1))).astype(float).tolist(),C2=(np.sum(c2[:,:,None]*u[k],axis=(0,1))/np.sum(abs(c2))).astype(float).tolist()))
 records.append(dict(n=n,rows=row,max_C1=max(abs(v) for q in row for v in q['C1']),max_C2=max(abs(v) for q in row for v in q['C2'])))
a.output.write_text(json.dumps(dict(scope=__doc__,records=records),indent=2)+'\n');print([(r['n'],r['max_C1'],r['max_C2']) for r in records])
