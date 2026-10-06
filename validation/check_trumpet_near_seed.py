"""Independent raw-tensor constraints at axes, horizon and trumpet interior."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'python'))
from hispid import Backend,Hole
from physical import constraints,norms

def run(library,only=None,factors=(.04,.02,.01)):
    b=Backend(str(library.resolve()));directions=np.r_[np.eye(3),[[.3,.7,.2]]];directions/=np.linalg.norm(directions,axis=1)[:,None]
    rows=[]
    for name,chi,v in [('spin99',.99,0.),('gamma10',0.,float(np.sqrt(.99)))]:
        if only and name!=only:continue
        h=Hole(1,spin=(0,0,chi),velocity=(v,0,0));G=1/np.sqrt(1-v*v);rh=np.sqrt(1-chi*chi)
        rest=np.concatenate([factor*rh*directions for factor in (.8,1.,1.5)])
        xyz=rest.copy();xyz[:,0]/=G
        step=np.linalg.norm(rest,axis=1)/G
        sample=lambda x:b.seed(h,x,choice=0,seed_family='trumpet_r0_m')
        row=dict(case=name,xyz=xyz.tolist(),constraints=[])
        for factor in factors:
            r=constraints(sample,xyz,factor*step)
            row['constraints'].append(dict(step_factor=factor,all=norms(r),interior=norms(r,np.arange(len(xyz))<4),horizon=norms(r,(np.arange(len(xyz))>=4)&(np.arange(len(xyz))<8)),exterior=norms(r,np.arange(len(xyz))>=8),minimum_metric_eigenvalue=float(np.min(r['min_metric_eigenvalue']))))
        row['passed']=bool(row['constraints'][-1]['all']['H_rms']<1e-7 and row['constraints'][-1]['all']['M_rms']<1e-7)
        rows.append(row)
    return dict(library_sha256=b.library_sha256(),rows=rows,passed=all(r['passed'] for r in rows),binary_acceptance=False,
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),ROOT/'validation/physical.py']})
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--library',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--case',choices=('spin99','gamma10'));ap.add_argument('--factors',default='.04,.02,.01');a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    factors=tuple(map(float,a.factors.split(',')))
    if not factors or any(not np.isfinite(h) or h<=0 for h in factors):raise ValueError('positive finite step factors required')
    r=run(a.library,a.case,factors);a.output.write_text(json.dumps(r,indent=2)+'\n')
    for row in r['rows']: print(row['case'],[(v['all']['H_rms'],v['all']['M_rms']) for v in row['constraints']],row['passed'])
