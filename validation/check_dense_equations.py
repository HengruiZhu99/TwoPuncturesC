"""Evaluate a retained polynomial on a denser grid, without another solve."""
import argparse,json
from pathlib import Path
import numpy as np
from hispid import Backend
from checkpoints import restore,select_record
from prolong import for_backend

p=argparse.ArgumentParser();p.add_argument('--library',required=True)
p.add_argument('--case',required=True);p.add_argument('--source-resolution',type=int,required=True)
p.add_argument('--source-nphi',type=int,required=True);p.add_argument('--dense-resolution',type=int,required=True)
p.add_argument('--dense-nphi',type=int,required=True);p.add_argument('--output',required=True)
p.add_argument('--dense-radial',type=int);p.add_argument('--dense-angular',type=int)
a=p.parse_args();backend=Backend(a.library)
rec=select_record(a.case,a.source_resolution,a.source_nphi);cfg,values=restore(backend,rec)
old=list(cfg.n);cfg.n[:]=[a.dense_radial or a.dense_resolution,a.dense_angular or a.dense_resolution,a.dense_nphi]
cfg.memory_limit_mib=8192;cfg.krylov_restart=32
values=for_backend(backend,values,old,list(cfg.n))
with backend.create(cfg) as s:
    s.set_unknowns(values);out=s.equation_samples()
g=out['attenuation'];psi=out['psi'];raw=out['physical_equivalent'].copy()
raw[:,0]*=-psi**5/8;raw[:,1:]*=psi[:,None]**10
groups={}
for label,mask in [('g0',g==0),('transition',(g>0)&(g<1)),('g1',g==1)]:
    groups[label]=dict(count=int(np.sum(mask)))
    if np.any(mask):
        groups[label].update(raw_conformal_linf=np.max(abs(raw[mask]),axis=0).tolist(),
            physical_equivalent_rms=np.sqrt(np.mean(out['physical_equivalent'][mask]**2,axis=0)).tolist(),
            physical_equivalent_linf=np.max(abs(out['physical_equivalent'][mask]),axis=0).tolist())
report=dict(case=a.case,source_resolution=old,evaluation_resolution=list(cfg.n),
    library_sha256=backend.library_sha256(),unknown_parameterization_id=backend.parameterization(),
    collocation_maps=backend.parameterization_maps(),source_diagnostics=rec['diagnostics'],groups=groups,
    solve_performed=False,acceptance=False,
    note='The same continuous retained modal polynomial is prolonged without solving. The denser native equations measure interpolation/PDE aliasing separately from an independent Cartesian finite-difference verifier. Modified interior equations do not establish physical vacuum.')
Path(a.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
