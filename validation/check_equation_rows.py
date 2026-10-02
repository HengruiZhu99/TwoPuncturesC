"""Inspect bound collocation residuals without loading another native context.

These are the modified interior equations, not vacuum physical constraints.
Recover their unweighted conformal values from the saved physical-equivalent
normalization and compare the inherited row weights by attenuation region.
"""
import argparse, json
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser();p.add_argument('--case',required=True)
p.add_argument('--output',required=True);a=p.parse_args()
root=Path(__file__).resolve().parents[1]
case=json.loads((root/'validation/results.json').read_text())[a.case]
records=[]
for rec in case['records']:
    n,nb,np_=rec['resolution']
    path=root/f"validation/raw/{rec['case']}_{n}_{np_}_collocation.npz"
    with np.load(path) as saved:
        xyz=saved['xyz'];g=saved['attenuation'];psi=saved['psi']
        raw=saved['physical_equivalent'].copy()
    raw[:,0]*=-psi**5/8;raw[:,1:]*=psi[:,None]**10
    angles=np.pi*(np.arange(n)+.5)/n
    polar=np.pi*(np.arange(nb)+.5)/nb
    w=np.broadcast_to((np.sin(polar)[:,None]*np.sin(angles)[None,:])**6,(np_,nb,n)).reshape(-1)
    weighted=w[:,None]*raw
    mass=sum(h['mass'] for h in rec['config']['hole'])
    distance=np.min([np.linalg.norm(xyz-np.array(h['center']),axis=1)/mass
                     for h in rec['config']['hole'] if h['mass']>0],axis=0)
    def stats(mask):
        if not np.any(mask):return dict(count=0)
        out=dict(count=int(np.sum(mask)),distance_range=[float(np.min(distance[mask])),float(np.max(distance[mask]))],
                 weight_range=[float(np.min(w[mask])),float(np.max(w[mask]))])
        for name,values in [('raw_conformal',raw),('weighted',weighted)]:
            values=values[mask]
            out[name]=dict(rms=np.sqrt(np.mean(values*values,axis=0)).tolist(),
                linf=np.max(abs(values),axis=0).tolist(),
                absolute_quantiles=np.quantile(abs(values),[.5,.9,.99,.999],axis=0).tolist())
        return out
    groups={label:stats(mask) for label,mask in [('g0',g==0),('transition',(g>0)&(g<1)),('g1',g==1)]}
    for lo,hi in [(0,.01),(.01,.03),(.03,.1),(.1,.3),(.3,1),(1,10),(10,np.inf)]:
        groups[f'distance_{lo}_{hi}M']=stats((distance>=lo)&(distance<hi))
    records.append(dict(resolution=rec['resolution'],library_sha256=rec['library_sha256'],
        unknown_parameterization_id=rec.get('unknown_parameterization_id'),
        collocation_maps=rec.get('collocation_maps'),diagnostics=rec['diagnostics'],groups=groups))
    print(json.dumps(dict(resolution=rec['resolution'],groups={k:groups[k] for k in ('g0','transition','g1')})),flush=True)
report=dict(case=a.case,records=records,acceptance=False,
    note='Diagnostic only. Interior equations are modified; no seed baseline is subtracted. Raw extrema and row-weight distributions do not replace independent exterior constraints or continuous horizon enclosure.')
Path(a.output).write_text(json.dumps(report,indent=2)+'\n')
