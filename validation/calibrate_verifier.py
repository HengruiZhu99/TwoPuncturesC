"""Measure the independent Cartesian verifier on exact vacuum controls.

Uses the saved binary's identical bulk points, then checks its three finest
grids at the same stencil sequence. This records a noise study; it does not
change or override the predeclared monotonic convergence gate.
"""
import argparse,json
import numpy as np
from hispid import Backend,Hole
from physical import constraints,norms
from run_validation import points
from checkpoints import ROOT,library_sha,restore

p=argparse.ArgumentParser();p.add_argument('--library',required=True)
p.add_argument('--case',default='highspin_stable')
p.add_argument('--steps',default='.016,.008,.004,.002,.001')
p.add_argument('--output',default='validation/verifier_floor_highspin.json')
args=p.parse_args();backend=Backend(args.library)
case=json.loads((ROOT/'validation/results.json').read_text())[args.case]
records=case['records'][-3:];cfg,_=restore(backend,records[-1])
x,near,bulk=points(cfg,records[-1]['horizon_scaled']);x=x[near:near+bulk]
steps=list(map(float,args.steps.split(',')))
output=dict(case=args.case,library_sha256=library_sha(backend),points=x.tolist(),
            steps=steps,controls=[],binary_grids=[],changes_original_gate=False)
def save():
    (ROOT/args.output).write_text(json.dumps(output,indent=2)+'\n')
def sequence(sample):
    values=[]
    for step in steps:
        result=constraints(sample,x,step)
        values.append(dict(step=step,norms=norms(result),H=result['H'].tolist(),M=result['M'].tolist()))
    return values
def brill_lindquist(x):
    x=np.asarray(x);psi=np.ones(len(x))
    for hole in cfg.hole:
        if hole.mass>0:psi+=hole.mass/(2*np.linalg.norm(x-np.array(hole.center),axis=1))
    return dict(gamma=(psi[:,None,None]**4*np.eye(3)).reshape(-1,9),Kij=np.zeros((len(x),9)),attenuation=np.ones(len(x)))
output['controls'].append(dict(case='exact analytic Brill-Lindquist',sequence=sequence(brill_lindquist)))
save();print('BL calibration saved',flush=True)
for i,hole in enumerate(cfg.hole):
    if hole.mass<=0:continue
    output['controls'].append(dict(case=f'exact Kerr seed {i}',
        sequence=sequence(lambda x,h=hole:backend.seed(h,x,cfg.conformal_choice))))
    save();print('seed calibration',i,'saved',flush=True)
for record in records:
    config,unknowns=restore(backend,record)
    with backend.create(config) as solution:
        solution.set_unknowns(unknowns)
        output['binary_grids'].append(dict(resolution=record['resolution'],sequence=sequence(solution.sample)))
    save();print('binary calibration',record['resolution'],'saved',flush=True)
print(json.dumps({k:[(r['step'],r['norms']['H_rms']) for r in v['sequence']]
                  for k,v in [(c['case'],c) for c in output['controls']]+[(str(c['resolution']),c) for c in output['binary_grids']]}))
