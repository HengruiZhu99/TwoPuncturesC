"""Independent physical FD at collocation nodes, to separate aliasing from mismatch.

No solve or native constraint formula is reused. Select a small fixed subset
nearest retained near/bulk off-grid observers and retain actual coordinates.
"""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'python'))
from hispid import Backend
from checkpoint_export import read_checkpoint
from physical import constraints,norms

def run(library,pilot,output):
    if output.exists():raise FileExistsError(output)
    b=Backend(str(library.resolve()));c,u,meta=read_checkpoint(pilot/'diagnostic.checkpoint')
    if meta['source_library_sha256']!=b.loaded_sha256 or meta['parameterization']!=b.parameterization():raise ValueError('checkpoint image/map mismatch')
    raw=np.load(pilot/'physical.npz');old=json.loads((pilot/'result.json').read_text());near=old['physical']['near']['count'];bulk=old['physical']['bulk']['count']
    selected=np.r_[np.argsort(raw['Mnorm'][:near])[-6:],near+np.argsort(raw['Mnorm'][near:near+bulk])[-6:]]
    A,B,P=map(int,c.n);maps=b.parameterization_maps();lam=maps['radial_stretch'];kap=maps['angular_stretch']
    sigma=.5*(1-np.cos(np.pi*(np.arange(A)+.5)/A));t=lam*sigma/(1-(1-lam)*sigma)
    eta=np.tanh(kap*-np.cos(np.pi*(np.arange(B)+.5)/B))/np.tanh(kap)
    centers=np.array([h.center for h in c.hole]);origin=centers.mean(axis=0);delta=centers[0]-centers[1];sep=np.linalg.norm(delta);axis=delta/sep
    e=np.eye(3)[np.argmin(abs(axis))];e=(e-axis*np.dot(e,axis));e/=np.linalg.norm(e);frame=np.array([axis,e,np.cross(axis,e)]).T
    axial=.5*sep*(1+t)[None,:]/(1-t)[None,:]*eta[:,None]
    rho=sep*np.sqrt(t)[None,:]/(1-t)[None,:]*np.sqrt(1-eta**2)[:,None]
    points=[];indices=[]
    for target in raw['xyz'][selected]:
        local=(target-origin)@frame;best=None
        for k in range(P):
            phi=2*np.pi*k/P;y=rho*np.cos(phi);z=rho*np.sin(phi)
            distance=(axial-local[0])**2+(y-local[1])**2+(z-local[2])**2;j,i=np.unravel_index(np.argmin(distance),distance.shape)
            if best is None or distance[j,i]<best[0]:best=(distance[j,i],(i,j,k),np.array([axial[j,i],y[j,i],z[j,i]])@frame.T+origin)
        indices.append([int(v) for v in best[1]]);points.append(best[2])
    points=np.array(points);distance=np.min(np.linalg.norm(points[:,None,:]-centers[None,:,:],axis=2),axis=1);step=np.minimum(.001,.001*distance)
    start=time.monotonic()
    with b.create_sampler(c) as s:
        s.set_unknowns(u);r=constraints(s.sample,points,step)
    np.savez_compressed(output.with_suffix('.npz'),xyz=points,step=step,**r)
    evidence=dict(max_unknown_by_mode_component=np.max(abs(u.reshape(P,B,A,4)),axis=(1,2)).tolist(),checkpoint=meta,driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),original_observer_indices=selected.tolist(),collocation_indices=indices,xyz=points.tolist(),step=step.tolist(),physical=dict(near=norms(r,np.arange(12)<6),bulk=norms(r,np.arange(12)>=6)),point_H=r['H'].tolist(),point_Mnorm=r['Mnorm'].tolist(),stencils_unmodified=bool(np.all(r['stencil_attenuation_all_one'])),seconds=time.monotonic()-start,binary_acceptance=False)
    output.write_text(json.dumps(evidence,indent=2)+'\n');print(json.dumps(evidence,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--pilot',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.library,a.pilot,a.output)
