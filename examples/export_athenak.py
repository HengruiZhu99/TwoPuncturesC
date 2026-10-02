"""Export a bound solved checkpoint or a fully specified exact seed control."""
import argparse,json
import numpy as np
from hispid import Backend,Hole
from checkpoint_export import write_checkpoint
from checkpoints import library_sha,select_record,restore

p=argparse.ArgumentParser();p.add_argument('--library',required=True);p.add_argument('--output',required=True)
group=p.add_mutually_exclusive_group(required=True)
group.add_argument('--case');group.add_argument('--seed',choices=('schwarzschild','kerr95','boost885','kerr95_boost885'))
p.add_argument('--resolution',type=int);p.add_argument('--nphi',type=int)
p.add_argument('--allow-diagnostic',action='store_true')
a=p.parse_args();b=Backend(a.library)
if a.seed:
    c=b.config();c.n[:]=[6,6,4];c.conformal_choice=0;c.inner_flatten=0
    c.omega[:]=[0,0];c.inner_min[:]=[0,0];c.inner_max[:]=[0,0];c.far_radius=0
    spin=(0,0,.95) if 'kerr95' in a.seed else (0,0,0)
    velocity=(.885,0,0) if 'boost885' in a.seed else (0,0,0)
    c.hole[0]=Hole(1,(0,0,0),spin,velocity);c.hole[1]=Hole(0,(-6,0,0))
    unknowns=np.zeros(4*np.prod(list(c.n)));acceptance='analytic_seed'
else:
    r=select_record(a.case,a.resolution,a.nphi);c,unknowns=restore(b,r)
    result=json.loads((__import__('checkpoints').ROOT/'validation/results.json').read_text())[a.case]
    converged=r['diagnostics']['status']==0
    if result.get('passed_strict') and r.get('passed_strict') and converged:acceptance='strong'
    elif result.get('passed') and r.get('passed_local') and converged:acceptance='preliminary'
    elif a.allow_diagnostic:acceptance='diagnostic'
    else:raise ValueError('case has not passed its physical convergence gate; --allow-diagnostic exports a labeled failure')
print(json.dumps(write_checkpoint(a.output,c,unknowns,library_sha(b),acceptance,b.parameterization()),indent=2))
