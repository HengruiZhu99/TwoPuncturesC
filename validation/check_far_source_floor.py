"""Calibrate isolated-vacuum source cancellation in a build's actual norm."""
import argparse,json,time
from pathlib import Path
import numpy as np
from hispid import Backend,Hole

def main():
    p=argparse.ArgumentParser();p.add_argument('--library',required=True);p.add_argument('--output',required=True)
    p.add_argument('--resolutions',default='32:32:16,64:64:28,80:160:28')
    a=p.parse_args();backend=Backend(a.library);records=[]
    spin=np.array([.2,-.3,.4]);spin*=.95/np.linalg.norm(spin)
    boost=np.array([-.7,.2,.3]);boost*=.885/np.linalg.norm(boost)
    output=dict(library_sha256=backend.library_sha256(),residual_scaling=backend.residual_scaling(),
        weighted_far_source_bound=1e-14,far_radius_minimum=100.,exact_correction=0.,records=records,passed=False,
        acceptance=False,note='Exact isolated vacuum seeds with all attenuation disabled; zero corrections. The far source floor is measured before any binary in a new norm. This does not accept a binary.')
    for shape in map(lambda x:list(map(int,x.split(':'))),a.resolutions.split(',')):
        for label,S,v in [('Schwarzschild',np.zeros(3),np.zeros(3)),('spin95',spin,np.zeros(3)),('boost885',np.zeros(3),boost),('combined_generic',spin,boost)]:
            cfg=backend.config();cfg.n[:]=shape;cfg.memory_limit_mib=8192;cfg.krylov_restart=32
            cfg.hole[0]=Hole(1,spin=S,velocity=v);cfg.hole[1]=Hole(0,center=(-6,0,0))
            cfg.conformal_choice=0;cfg.inner_flatten=0;cfg.inner_min[:]=cfg.inner_max[:]=[0,0]
            cfg.omega[:]=[0,0];cfg.far_radius=0;start=time.monotonic()
            with backend.create(cfg) as solution:
                weighted=solution.residual(solution.unknowns()).reshape(-1,4)
                sampled=solution.equation_samples()
            radius=np.linalg.norm(sampled['xyz'],axis=1);mask=radius>=100
            if not np.any(mask):raise ValueError('no far nodes')
            maxima=np.max(abs(weighted[mask]),axis=0)
            row=dict(case=label,resolution=shape,far_node_count=int(mask.sum()),far_radius_range=[float(radius[mask].min()),float(radius[mask].max())],
                weighted_far_source_linf=maxima.tolist(),physical_equivalent_far_linf=np.max(abs(sampled['physical_equivalent'][mask]),axis=0).tolist(),
                passed=bool(np.all(np.isfinite(maxima)) and maxima.max()<1e-14),seconds=time.monotonic()-start)
            records.append(row);Path(a.output).write_text(json.dumps(output,indent=2)+'\n');print(label,shape,maxima,row['passed'],flush=True)
    output['passed']=all(r['passed'] for r in records);Path(a.output).write_text(json.dumps(output,indent=2)+'\n')
    if not output['passed']:raise SystemExit(1)

if __name__=='__main__':main()
