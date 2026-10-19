"""Fixed-equation low-order studies; all failures remain diagnostic results.

The paper binary has m+=m-=b=1, P+=-P-=0.2 along y (Fig.4).
The Kerr analogue matches isolated charges, not the BY physical free metric.
No comparison between the two formulations asserts identical physical data.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'python'), str(ROOT/'examples'), str(ROOT/'validation')]
from hispid import Backend as HiSpID, Hole
from bowen_york import Backend as BowenYork
from configs import as_dict
from trumpet_configs import qi_moderate_wide
from physical import constraints, norms
from run_validation import points
from parallel_observer import ParallelObserver
from prolong import for_backend
from checkpoint_export import write_checkpoint


def configuration(backend, case, n):
    if case == 'moderate':
        return qi_moderate_wide(backend, n, 16)
    c = backend.config()
    c.n[:] = [n, n, max(4, 2*((n+3)//4))]
    c.conformal_choice = 0
    c.inner_flatten = 1
    c.far_radius = 0
    c.omega[:] = [0.5, 0.5]
    v = .2/np.sqrt(1+.2**2)
    for i, sign in enumerate((1, -1)):
        c.hole[i] = Hole(1, (sign, 0, 0), velocity=(0, sign*v, 0))
        c.inner_min[i] = .075/np.sqrt(1+.2**2)
        c.inner_max[i] = .425/np.sqrt(1+.2**2)
    c.tolerance = 1e-12
    c.max_newton = 24
    return c


def by_config(backend, c):
    q = backend.config()
    q['integer'].update(zip(('npoints_A', 'npoints_B', 'npoints_phi'), list(c.n)))
    q['integer'].update(Newton_maxit=24, verbose=1)
    q['real'].update(par_b=float(c.hole[0].center[0]), Newton_tol=c.tolerance)
    for side, h in zip(('plus', 'minus'), c.hole):
        v, S = np.array(h.velocity), np.array(h.spin)
        G = 1/np.sqrt(1-v@v)
        P = h.mass*G*v
        S = G*S-G*G/(G+1)*(v@S)*v
        q['real']['par_m_'+side] = h.mass
        for k in range(3):
            q['real']['par_P_'+side+str(k+1)] = float(P[k])
            q['real']['par_S_'+side+str(k+1)] = float(S[k])
    return q


def differences(a, b):
    out = {}
    for name in ('gamma', 'Kij', 'psi'):
        delta = np.asarray(a[name])-np.asarray(b[name])
        out[name] = dict(absolute_linf=float(np.max(abs(delta))),
                         scaled_linf=float(np.max(abs(delta)/(1+abs(b[name])))))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hispid-library', type=Path, required=True)
    p.add_argument('--by-library', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--case', choices=('moderate', 'paper-binary'), required=True)
    p.add_argument('--levels', default='32,48,64')
    p.add_argument('--execution', choices=('reference', 'kokkos'), default='kokkos')
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--reference-hispid', type=Path)
    p.add_argument('--reference-by', type=Path)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    levels = list(map(int, a.levels.split(',')))
    if levels != sorted(set(levels)) or min(levels) < 8:
        raise ValueError('strictly increasing grid counts >=8 required')
    hb, bb = HiSpID(str(a.hispid_library.resolve())), BowenYork(str(a.by_library.resolve()))
    c0 = configuration(hb, a.case, levels[0])
    xyz, near, bulk = points(c0, True)
    xyz = xyz[:near+bulk]
    step = 4*np.minimum(.002, .002*np.min([np.linalg.norm(xyz-np.array(h.center), axis=1) for h in c0.hole], axis=0))
    report = dict(case=a.case, levels=levels, records=[],
                  reference='https://arxiv.org/pdf/gr-qc/0404056',
                  scope='same solution within each fixed formulation; BY and Kerr free metrics differ',
                  observer=dict(order=6, xyz=xyz.tolist(), step=step.tolist()),
                  libraries=dict(hispid=hb.library_sha256(), bowen_york=bb.library_sha256()),
                  hispid_basis=dict(parameterization=hb.parameterization(), maps=hb.parameterization_maps()),
                  driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    previous = {}
    saved = {}
    prev_values = prev_shape = None
    def save():
        (a.output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    save()
    for n in levels:
        c = configuration(hb, a.case, n)
        c.memory_limit_mib = 65536
        c.krylov_restart = 64
        c.max_krylov = 512
        for system in ('bowen_york', 'hispid'):
            directory = a.output/f'{system}{n}'
            directory.mkdir()
            start = time.monotonic()
            rec = dict(system=system, n=list(c.n), config=as_dict(c))
            if system == 'bowen_york':
                solution = bb.create(by_config(bb, c))
            else:
                solution = hb.create(c, execution=a.execution, geometry='host')
            with solution as s:
                if system == 'hispid' and prev_values is not None:
                    s.set_unknowns(for_backend(hb, prev_values, prev_shape, list(c.n)))
                rec['diagnostics'] = (s.solve(linear_rtol=.001, krylov='gmres', preconditioner='modal',
                    max_krylov=512, restart=64, execution=a.execution, memory_limit_mib=65536)
                    if system == 'bowen_york' else s.solve(linear_rtol=.1, krylov='gmres'))
                rec['solve_seconds'] = time.monotonic()-start
                rec['work_statistics'] = s.work_statistics()
                values = s.unknowns()
                np.save(directory/'unknowns.npy', values)
                if system == 'hispid':
                    write_checkpoint(directory/'diagnostic.checkpoint', c, values, hb.library_sha256(), 'diagnostic', hb.parameterization())
                    prev_values, prev_shape = values, list(c.n)
                    with ParallelObserver(s, a.threads) as observer:
                        fields = observer(xyz)
                        physical = constraints(observer, xyz, step, order=6)
                else:
                    fields = s.sample(xyz)
                    physical = constraints(s.sample, xyz, step, order=6)
                rec['physical'] = {name: norms(physical, mask) for name, mask in
                    [('near', np.arange(len(xyz)) < near), ('bulk', np.arange(len(xyz)) >= near)]}
                np.savez_compressed(directory/'fields.npz', xyz=xyz, psi=fields['psi'], gamma=fields['gamma'], Kij=fields['Kij'], H=physical['H'], M=physical['M'])
                saved[system, n] = {k:fields[k].copy() for k in ('psi', 'gamma', 'Kij')}
                if system in previous:
                    rec['difference_from_previous'] = differences(fields, previous[system])
                previous[system] = saved[system, n]
                refpath = a.reference_hispid if system == 'hispid' else a.reference_by
                if refpath:
                    with np.load(refpath) as z:
                        if not np.array_equal(z['xyz'][:len(xyz)], xyz):
                            raise ValueError('reference observation points differ')
                        rec['difference_from_existing_fine'] = {
                            k:dict(absolute_linf=float(np.max(abs(fields[k]-z[k][:len(xyz)].reshape(fields[k].shape)))),
                                   scaled_linf=float(np.max(abs(fields[k]-z[k][:len(xyz)].reshape(fields[k].shape))/(1+abs(z[k][:len(xyz)].reshape(fields[k].shape)))))) for k in ('gamma', 'Kij')}
                if system == 'hispid' and n == levels[-1]:
                    # Independent initial guess, unchanged equations/grid. Do
                    # not mistake agreement between different formulations.
                    with hb.create(c, execution=a.execution, geometry='host') as cold:
                        cold_diag = cold.solve(linear_rtol=.1, krylov='gmres')
                        rec['cold_start'] = dict(diagnostics=cold_diag, difference=differences(cold.sample(xyz), fields))
            if system == 'bowen_york' and n == levels[0]:
                with bb.create(by_config(bb, c)) as baseline:
                    diag = baseline.solve(preconditioner='lines', execution='reference', max_krylov=100)
                    rec['vanilla_control'] = dict(diagnostics=diag, difference=differences(baseline.sample(xyz), fields))
            rec['total_seconds'] = time.monotonic()-start
            report['records'].append(rec)
            save()
            print(system, n, rec['diagnostics'], rec['physical'], flush=True)
    for rec in report['records']:
        system = rec['system']
        if rec['n'][0] != levels[-1]:
            rec['difference_from_study_finest'] = differences(saved[system, rec['n'][0]], saved[system, levels[-1]])
    save()


if __name__ == '__main__':
    main()
