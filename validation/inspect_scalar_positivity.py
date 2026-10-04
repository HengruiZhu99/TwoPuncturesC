"""Fixed-state scalar positivity diagnosis; no native image or solve is loaded.

Scope: unboosted QI seeds, symmetric x-axis foci, current mapped modal basis,
and no far filter. Barycentric P evaluation bypasses native coefficient sums.
The QI scalar is evaluated from the analytic Kerr radial formula. This is a
limited reconstruction diagnosis, never independent vacuum acceptance.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from prolong import chebyshev_evaluation_matrix


def seed_scalar(x,config):
    result=np.ones(len(x))
    for hole in config['hole']:
        mass=hole['mass'];spin=np.asarray(hole['spin']);d=x-np.asarray(hole['center'])
        r=np.linalg.norm(d,axis=1);a=np.linalg.norm(spin)/mass
        if np.any(r==0):raise ValueError('puncture excluded from scalar diagnosis')
        costheta=np.sum(d*(spin/np.linalg.norm(spin))[None,:],axis=1)/r if a else np.zeros(len(x))
        rb=r+mass+(mass*mass-a*a)/(4*r)
        result+=((rb*rb+a*a*costheta*costheta)/(r*r))**.25-1
    return result


def scalar_correction(x,values,shape,maps,b):
    na,nb,nphi=shape;lam,kappa=maps['radial_stretch'],maps['angular_stretch']
    radii=[np.linalg.norm(x-np.array([sign*b,0,0]),axis=1) for sign in (-1,1)]
    total=sum(radii);t=(total-2*b)/(total+2*b);eta=2*x[:,0]/total
    sigma=t/(lam+(1-lam)*t);zeta=np.arctanh(eta*np.tanh(kappa))/kappa
    ea=chebyshev_evaluation_matrix(na,2*sigma-1)
    eb=chebyshev_evaluation_matrix(nb,zeta)
    grid=values.reshape(nphi,nb,na,4)[:,:,:,0]
    mode=np.einsum('qi,mji->qmj',ea,grid,optimize=True)
    mode=np.einsum('qj,qmj->qm',eb,mode,optimize=True)
    phi=np.arctan2(x[:,2],x[:,1]);sn=np.sqrt(np.maximum(0,1-eta*eta));a=np.sqrt(t)
    result=np.zeros(len(x));half=nphi//2
    for k in range(nphi):
        m=k if k<=half else k-half;r=m if m<=4 else (3 if m%2 else 4)
        norm=np.sqrt((1 if m in (0,half) else 2)/nphi)
        phase=np.cos(m*phi) if k<=half else np.sin(m*phi)
        result+=-2*(1-t)*(a*sn)**r*norm*phase*mode[:,k]
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('result','case','solve-artifact','output'):p.add_argument('--'+name,required=True)
    a=p.parse_args();parent=Path(a.result).resolve(strict=True);raw=Path(a.solve_artifact).resolve(strict=True)
    out=Path(a.output).resolve()
    if out.exists():raise FileExistsError('preserve earlier scalar diagnosis')
    sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    group=json.loads(parent.read_text())[a.case];row=group['incomplete_attempt']['record'];cfg=row['config']
    if (row['unknown_parameterization_id']!='modal_P_C2prolate_mapped_v2'
        or cfg['conformal_choice']!=0 or cfg['far_radius']!=0
        or any(np.any(h['velocity']) or h['mass']<=0 for h in cfg['hole'])
        or cfg['hole'][0]['center'][1:]!=[0.,0.] or cfg['hole'][1]['center']!=[-cfg['hole'][0]['center'][0],0.,0.]):
        raise ValueError('case outside the declared limited reconstruction scope')
    bound={str(parent):sha(parent),str(raw):sha(raw),str(Path(__file__).resolve()):sha(__file__),
        str(Path(__import__('prolong').__file__).resolve()):sha(__import__('prolong').__file__)}
    if bound[str(raw)]!=row['solve_artifact']['sha256']:raise ValueError('retained solve coefficients changed')
    values=np.load(raw,allow_pickle=False)['unknowns'];shape=cfg['n'];maps=row['collocation_maps'];b=cfg['hole'][0]['center'][0]
    if values.shape!=(4*np.prod(shape),) or not np.isfinite(values).all():raise ValueError('invalid retained unknowns')
    dirs=np.array([[.73,.31,.61],[-.41,.82,.39],[.22,-.51,.83],[-.69,-.44,.57],[.39,.73,-.56],[.81,-.38,-.45]])
    dirs/=np.linalg.norm(dirs,axis=1)[:,None];near=[]
    for hole in cfg['hole']:
        mass=hole['mass'];chi=np.linalg.norm(hole['spin'])/mass**2;rh=.5*mass*np.sqrt(1-chi*chi)
        for r in (1.5*rh,3*rh,.6*mass,mass):near.extend(np.asarray(hole['center'])+r*dirs)
    x=np.r_[np.array(near),dirs*4,dirs*8,dirs*16]
    steps=np.asarray(row['verifier_steps'])
    if len(x)!=len(steps):raise ValueError('retained point/step inventory differs')
    offsets=[np.zeros(3)]
    for d in range(3):
        for w in (-2,-1,1,2):v=np.zeros(3);v[d]=w;offsets.append(v)
    for d in range(3):
        for e in range(d+1,3):
            for u in (-2,-1,1,2):
                for v in (-2,-1,1,2):q=np.zeros(3);q[d]=u;q[e]=v;offsets.append(q)
    offsets=np.asarray(offsets);points=(x[None,:,:]+offsets[:,None,:]*steps[None,:,None]).reshape(-1,3)
    u=np.concatenate([scalar_correction(block,values,shape,maps,b) for block in np.array_split(points,32)])
    seed=seed_scalar(points,cfg);psi=seed+u;worst=np.argsort(psi)[:8]
    if not all(np.isfinite(v).all() for v in (points,u,seed,psi)):
        raise ValueError('nonfinite scalar reconstruction; no positivity conclusion')
    result=dict(purpose='fixed_iterate_scalar_positivity_diagnosis',solve_performed=False,physical_acceptance=False,
        source_library_sha256=row['library_sha256'],config=cfg,maps=maps,bound_artifacts_sha256=bound,
        native_image_loaded=False,all_values_finite=True,method='barycentric mapped modal P and analytic unboosted Kerr QI scalar',
        point_count=len(x),stencil_point_count=len(points),negative_or_zero_count=int(np.sum(psi<=0)),
        min_scalar=float(np.min(psi)),center_sample_min_scalar=float(np.min(psi[:len(x)])),
        worst=[dict(index=int(i),sample_index=int(i%len(x)),offset=offsets[i//len(x)].tolist(),
            xyz=points[i].tolist(),seed_scalar=float(seed[i]),correction=float(u[i]),scalar=float(psi[i])) for i in worst])
    if any(sha(path)!=s for path,s in bound.items()):raise ValueError('bound positivity diagnosis inputs changed')
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('negative_or_zero_count','min_scalar','center_sample_min_scalar','worst')}))


if __name__=='__main__':main()
