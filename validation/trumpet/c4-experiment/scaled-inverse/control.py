"""Independent-process equation witness and optional small nonlinear control."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(ROOT / 'python'), str(ROOT / 'examples')]
from hispid import Backend
from trumpet_configs import trumpet_moderate

p = argparse.ArgumentParser()
p.add_argument('--library', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--solve', action='store_true')
p.add_argument('--krylov', choices=['gmres','bicgstab','lgmres'], default='gmres')
p.add_argument('--krylov-restart', type=int)
p.add_argument('--execution', choices=['reference','kokkos'], default='reference')
a = p.parse_args()
if a.output.exists() or a.output.with_suffix('.npz').exists():
    raise FileExistsError(a.output)
b = Backend(str(a.library.resolve()))
c = trumpet_moderate(b, 12, 16)
c.memory_limit_mib = 2048
if a.krylov_restart is not None:c.krylov_restart=a.krylov_restart
with b.create(c, execution=a.execution, geometry='host') as s:
    k = np.arange(s.size)
    base, direction = 1e-8*np.sin(.13*k), 1e-7*np.cos(.27*k)
    residual, jvp = s.residual(base), s.jvp(base, direction)
    assert np.isfinite(residual).all() and np.isfinite(jvp).all()
    np.savez_compressed(a.output.with_suffix('.npz'), residual=residual, jvp=jvp)
    result = dict(execution=a.execution, library_sha256=b.loaded_sha256, parameterization=b.parameterization(),
                  residual_scaling=b.residual_scaling(), binary_acceptance=False,
                  driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if a.solve:
        result['diagnostics'] = s.solve(krylov=a.krylov, linear_rtol=.001)
        result['resolved_options'] = s.resolved_options
        result['linear_history'] = s.linear_history()
        if result['diagnostics']['converged']:
            xyz = np.array([[3.3,.4,.2],[2.7,-.5,.3],[-3.2,.3,-.4],[-2.7,-.4,-.3],[.3,1.2,.7],[4,2,-3]])
            result['field_witness'] = {k:v.tolist() for k,v in s.sample_with_derivatives(xyz).items()}
a.output.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
