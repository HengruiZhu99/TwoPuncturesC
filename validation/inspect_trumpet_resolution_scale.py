"""Count meridional nodes resolving the fixed interior transition shells."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np

def inspect(path,radial_stretch,angular_stretch):
    d=json.loads(path.read_text());c=d['config'];na,nb,nphi=c['n']
    centers=np.array([h['center'] for h in c['hole']])
    if not np.all(centers[:,1:]==0) or not centers[0,0]>centers[1,0]:
        raise ValueError('this diagnostic expects holes ordered along the x axis')
    b=(centers[0,0]-centers[1,0])/2
    sigma=.5*(1-np.cos(np.pi*(np.arange(na)+.5)/na))
    t=radial_stretch*sigma/(1-(1-radial_stretch)*sigma)
    eta=np.tanh(angular_stretch*(-np.cos(np.pi*(np.arange(nb)+.5)/nb)))/np.tanh(angular_stretch)
    # Exact prolate distance to the +x focus: r+=b(cosh(X)-cos(R)).
    radius=b*((1+t[:,None])/(1-t[:,None])-eta[None,:])
    radial_axis=2*b*t/(1-t);polar_axis=b*(1-eta)
    holes=[]
    for h in range(2):
        lo,hi=c['inner_min'][h],c['inner_max'][h];r=radius if h==0 else radius[:,::-1]
        count=lambda a:int(np.count_nonzero((a>lo)&(a<hi)))
        holes.append(dict(hole=h,inner_min=lo,inner_max=hi,
            radial_axis_nodes_in_transition=count(radial_axis),polar_axis_nodes_in_transition=count(polar_axis),
            meridional_nodes_in_transition=count(r)))
    return dict(run=str(path.parent),input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),resolution=c['n'],
        radial_stretch=radial_stretch,angular_stretch=angular_stretch,holes=holes)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',type=Path,nargs='+',required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    rows=[inspect(run/'result.json',r,k) for run in a.runs for r,k in ((.2,2.),(.03,3.5))]
    a.output.write_text(json.dumps(dict(rows=rows,driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        note='Geometric node counts only. Axis projections are not actual endpoint collocation nodes. The candidate map is not a solved or selected configuration.'),indent=2)+'\n')
