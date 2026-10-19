"""One finest-grid perturbed-start and independent FD observer check."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Backend
from checkpoint_export import read_checkpoint
from parallel_observer import ParallelObserver
from physical import constraints,norms

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',type=Path,required=True)
p.add_argument('--checkpoint',type=Path,required=True);p.add_argument('--study',type=Path,required=True)
p.add_argument('--output',type=Path,required=True);p.add_argument('--constant-p-start',action='store_true',help='diagnose an initial guess that violates homogeneous focus rows');a=p.parse_args()
b=Backend(str(a.library.resolve()));c,u,metadata=read_checkpoint(a.checkpoint)
study=json.loads(a.study.read_text());x=np.array(study['observer']['xyz']);step=np.array(study['observer']['steps'])
with b.create_sampler(c) as reference,b.create(c,execution='kokkos') as s:
    reference.set_unknowns(u);original=reference.sample(x)
    du=np.tile([1e-5,1e-6,-2e-6,3e-6],len(u)//4) if a.constant_p_start else .05*u;s.set_unknowns(u+du)
    diagnostics=s.solve(krylov='gmres');final=s.sample(x)
    difference={k:float(np.max(abs(final[k]-original[k]))) for k in ('gamma','Kij')}
    with ParallelObserver(s,16) as sample:
        physical=constraints(sample,x,step,6);half=constraints(sample,x,step/2,6)
    observer_change={k:float(np.max(abs(physical[k]-half[k]))) for k in ('H','M')}
    result=dict(checkpoint=metadata,library_sha256=b.library_sha256(),diagnostics=diagnostics,
        start_perturbation=([1e-5,1e-6,-2e-6,3e-6] if a.constant_p_start else 'scale the saved unknowns by1.05; preserve homogeneous focus/axis conditions'),final_field_difference=difference,
        final_unknown_linf=float(np.max(abs(s.unknowns()-u))),
        independent_fd6=norms(physical),half_step_fd6=norms(half),observer_change=observer_change,
        passed=bool(diagnostics['converged'] and max(difference.values())<1e-9 and max(observer_change.values())<1e-8))
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    assert result['passed'],result
