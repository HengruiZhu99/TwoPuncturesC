"""Recheck retained numerical evidence without repeating or replacing solves."""
import argparse
import json
import hashlib
import subprocess
from pathlib import Path
import numpy as np
from benchmark_bowen_york import ROOT,digest


def compare_states(first, second, shared=False):
    with np.load(first) as a,np.load(second) as b:
        if not shared and set(a.files)!=set(b.files):raise AssertionError('array keys differ')
        checks={}
        for name in a.files:
            x,y=a[name],b[name]
            checks[name]=dict(shape=list(x.shape),dtype=str(x.dtype),
                              bitwise_equal=x.shape==y.shape and x.dtype==y.dtype and x.tobytes()==y.tobytes(),
                              finite=bool(np.isfinite(x).all() and np.isfinite(y).all()))
        if not all(c['bitwise_equal'] and c['finite'] for c in checks.values()):raise AssertionError(checks)
        return checks


def verify_evidence(evidence):
    checks={}
    for mode,versions in evidence['records'].items():
        reference=versions['before'][0];checks[mode]=[]
        trace=lambda r:[line for line in (ROOT/r['log']).read_text().splitlines() if line.startswith(('Newton:','bicgstab:'))]
        for version,series in versions.items():
            for r in series:
                parameters_equal=r['config']==reference['config'] and all(r.get(k)==reference.get(k) for k in ('by_real_parameters','by_integer_parameters'))
                diagnostic=lambda record:{k:v for k,v in record['diagnostics'].items() if k!='seconds'}
                diagnostics_equal=diagnostic(r)==diagnostic(reference)
                traces_equal=bool(trace(r)) and trace(r)==trace(reference) if mode=='by' else True
                fingerprints_current=digest(ROOT/r['state'])==r['state_sha256'] and digest(ROOT/r['log'])==r['log_sha256']
                if not (parameters_equal and diagnostics_equal and traces_equal and fingerprints_current and r['converged']):
                    raise RuntimeError('configuration, diagnostics, trace, raw hash or native convergence verification failed')
                arrays=compare_states(ROOT/reference['state'],ROOT/r['state'])
                checks[mode].append(dict(version=version,parameters_equal=parameters_equal,diagnostics_equal=diagnostics_equal,
                                         trace_equal=traces_equal,raw_fingerprints_current=fingerprints_current,arrays=arrays))
    return checks


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',required=True);p.add_argument('--output',required=True)
    args=p.parse_args();out=Path(args.output)
    if out.exists():raise FileExistsError(out)
    evidence=json.loads(Path(args.evidence).read_text());checks=verify_evidence(evidence)
    raw=ROOT/'validation/raw/by_upstream_controls';original={}
    upstream40=json.loads((raw/'upstream40.json').read_text())
    if not upstream40['converged'] or digest(raw/'upstream40.npz')!=upstream40['state_sha256']:
        raise RuntimeError('original full-grid convergence/state fingerprint failed')
    original['moderate_40']=dict(upstream_commit='ec563aeb672235b9443c330f9cde65f7246e8ea4',
                               arrays=compare_states(raw/'upstream40.npz',ROOT/evidence['records']['by']['after'][0]['state'],shared=True),
                               converged=json.loads((raw/'upstream40.json').read_text())['converged'])
    for case in ('moderate_small','target_mass','bl','spin95'):
        a=json.loads((raw/(case+'_upstream.json')).read_text());b=json.loads((raw/(case+'_optimized.json')).read_text())
        if not (a['converged'] and b['converged']):raise RuntimeError('original small-grid native convergence failed')
        diagnostic_keys=('residual_linf','final_bare_masses','internal_end_adm_masses','target_mass_errors','real_parameters','integer_parameters')
        if not all(a[k]==b[k] for k in diagnostic_keys):raise RuntimeError('original small-grid configuration/diagnostics differ')
        for version,r in (('upstream',a),('optimized',b)):
            if digest(raw/(case+'_'+version+'.npz'))!=r['state_sha256']:raise RuntimeError('original small-grid state fingerprint failed')
        original[case]=dict(grid=[12,18,8],arrays=compare_states(raw/(case+'_upstream.npz'),raw/(case+'_optimized.npz')),
                            upstream=a,optimized=b,diagnostics_equal=True)
    files=['src/TP_Newton.c','include/TP_LineCache.h','src/HiSpID_solver.cpp','src/HiSpID_geometry.cpp',
           'tests/test_by_line_cache.c','tests/test_hispid_krylov_memory.cpp','validation/benchmark_bowen_york.py',
           'validation/benchmark_solver_efficiency.py','validation/compare_by_upstream.py','validation/verify_solver_efficiency.py',
           'build-hispid-efficiency/baseline/newton.c','build-hispid-efficiency/baseline/timer.c',
           'build-hispid-efficiency/optimized/newton.c','build-hispid-efficiency/optimized/timer.c',
           'build-hispid-efficiency/final-tests.log','build-hispid-efficiency/final-physical-api.log']
    files += [str(f.relative_to(ROOT)) for f in raw.iterdir() if f.is_file()]
    libraries={}
    for mode,versions in evidence['records'].items():
        for version,series in versions.items():
            # Resolve explicit CLI library path from the preserved command.
            flag='--hispid-library' if mode=='hispid' else '--by-library';command=series[0]['command'];path=command[command.index(flag)+1]
            if digest(path)!=series[0]['library_sha256']:raise RuntimeError('benchmark library fingerprint changed')
            libraries[mode+'_'+version]=dict(path=path,sha256=digest(path))
    original_root=ROOT/'build-hispid-efficiency/upstream/source'
    original_sources={}
    for relative in ('src/TP_Newton.c','src/TP_Equations.c','src/TP_FuncAndJacobian.c',
                     'src/TP_CoordTransf.c','src/TP_Utilities.c','src/TwoPunctures.c','src/TwoPuncturesRun.c','include/TwoPunctures.h'):
        path=original_root/relative
        committed=subprocess.check_output(['git','show','ec563aeb:'+relative],cwd=ROOT)
        expected=hashlib.sha256(committed).hexdigest()
        if digest(path)!=expected:raise RuntimeError('original source differs from ec563aeb: '+relative)
        original_sources[relative]=expected
    upstream_library=original_root/'lib/libTwoPunctures.so'
    if digest(upstream_library)!=upstream40['library_sha256']:raise RuntimeError('original full-grid library fingerprint changed')
    libraries['original_standalone']=dict(path=str(upstream_library),sha256=digest(upstream_library),source_sha256=original_sources)
    final_library=ROOT/'build-hispid-efficiency-final/lib/libTwoPunctures.so'
    if digest(final_library)!=original['target_mass']['optimized']['library_sha256']:raise RuntimeError('final BY library fingerprint changed')
    libraries['optimized_clean_build']=dict(path=str(final_library),sha256=digest(final_library))
    result=dict(evidence=args.evidence,evidence_sha256=digest(args.evidence),passed=True,
                strengthened_retained_record_checks=checks,original_standalone_controls=original,
                provenance=dict(source_and_raw_sha256={f:digest(ROOT/f) for f in files},libraries=libraries,
                                baseline_commit='c4158cb2be0fcffb7db43b42c19c7ada1d556214',
                                before_hispid_library_source='9f85d08264fd04ca30a4d61201c35302d2892765',
                                before_by_core_matches_upstream='ec563aeb672235b9443c330f9cde65f7246e8ea4',
                                compiler='Apple clang21.0.0',GSL='2.8',flags='C99/C++17 -O3 -fPIC',
                                trace_generation='Before: git show c4158cb:src/TP_Newton.c; after: current TP_Newton.c. Replace %10.3e with %.17e and |F|=%e with |F|=%.17e for identical full-precision printf-only traces. Compile by_newton_timer.c pointing at each snapshot instead of TP_Newton.o; link each matching C object bank.'),
                note='This verification reads retained records; it does not overwrite or rerun historical measurements. Source/script hashes describe the final replay tools; the compiled native image fingerprints and raw snapshot/trace hashes identify the actual numerical witnesses. Current target-mass fixtures test the first fresh solve, not original static-lifecycle reuse. Spin95 here labels a BY input control, not HiSpID binary physical validation.')
    out.write_text(json.dumps(result,indent=2)+'\n');print('retained-record and original standalone verification passed',flush=True)
if __name__=='__main__':main()
