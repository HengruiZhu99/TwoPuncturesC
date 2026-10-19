"""Same single-hole BY data on spherical and unchanged binary native maps."""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'python'), str(ROOT/'validation')]
from bowen_york import Backend
from single_hole import BowenYorkSingleHole

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--library', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--execution', choices=('reference', 'kokkos'), default='kokkos')
a = p.parse_args()
if a.output.exists():
    raise FileExistsError(a.output)
b = Backend(str(a.library.resolve()))
xyz = np.array([[.5, .1, .2], [1., .2, -.4], [3., 2., 4.], [.17, -.12, .1], [8., -3., 2.]])
report = dict(scope=__doc__, library_sha256=b.library_sha256(), records=[])
for label, spin, momentum in (('spin', .2, 0.), ('boost', 0., .2)):
    direct = b.single_hole(spin=spin, momentum=momentum, nr=100, ntheta=32)
    direct_diagnostics = direct.solve()
    truth = direct.sample(xyz)
    for n in (32, 64):
        c = b.config()
        c['integer'].update(npoints_A=n, npoints_B=2*n, npoints_phi=16, Newton_maxit=16)
        c['real'].update(par_m_plus=1., par_m_minus=0., center_offset1=-3.,
                         par_S_plus3=spin, par_P_plus3=momentum)
        with b.create(c) as s:
            diag = s.solve(linear_rtol=.001, krylov='gmres', preconditioner='modal', max_krylov=512, restart=64, execution=a.execution)
            actual = s.sample(xyz)
            difference = {k:float(np.max(abs(actual[k]-truth[k])/(1+abs(truth[k])))) for k in ('gamma', 'Kij', 'psi')}
        report['records'].append(dict(case=label, n=[n, 2*n, 16], native=diag, direct=direct_diagnostics,
                                      difference=difference, config=c))
        a.output.write_text(json.dumps(report, indent=2)+'\n')
        print(label, n, diag, difference, flush=True)
report['passed'] = bool(all(r['direct']['converged'] and max(r['difference'].values()) < 2e-7
    for r in report['records'] if r['n'][0] == 64))
report['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
a.output.write_text(json.dumps(report, indent=2)+'\n')
if not report['passed']:
    raise RuntimeError('native and direct single-hole data differ')
