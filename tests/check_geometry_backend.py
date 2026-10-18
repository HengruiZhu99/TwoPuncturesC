"""One small public-API parity check for the selected production geometry path."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'python'),str(ROOT/'examples')]
from hispid import Backend
from trumpet_configs import qi_moderate_wide
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();b=Backend(str(a.library.resolve()));c=qi_moderate_wide(b,12,16);c.memory_limit_mib=2048;c.krylov_restart=16
n=4*np.prod(c.n);x=1e-8*np.cos(.23*np.arange(n));v=1e-7*np.sin(.31*np.arange(n))
with b.create(c) as ref:
 r0=ref.residual(x);j0=ref.jvp(x,v)
with b.create(c,execution='kokkos',geometry='host') as dev:
 r1=dev.residual(x);j1=dev.jvp(x,v)
relative=lambda a,b:float(np.linalg.norm(a-b)/(1e-30+np.linalg.norm(a)))
result=dict(residual_relative_l2=relative(r0,r1),jvp_relative_l2=relative(j0,j1),library_sha256=b.library_sha256(),parameterization=b.parameterization())
assert max(result['residual_relative_l2'],result['jvp_relative_l2'])<1e-10,result
result['passed']=True;a.output.write_text(json.dumps(result,indent=2)+'\n');print(result)
