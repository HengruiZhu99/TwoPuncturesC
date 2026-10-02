"""Independent retained-data controls plus original standalone small fixtures.

Requires accepted common-stopping measurements; new numerical workers run serially.
"""
import argparse,json,os,subprocess,sys
from pathlib import Path
import numpy as np
from benchmark_bowen_york import ROOT,digest
from verify_solver_efficiency import compare_states


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',required=True);p.add_argument('--by-library',required=True);p.add_argument('--output',required=True);args=p.parse_args()
    out=Path(args.output);raw=ROOT/'validation/raw'/out.stem
    if out.exists():raise FileExistsError(out)
    evidence=json.loads(Path(args.evidence).read_text())
    if not evidence['passed']:raise RuntimeError('shared stopping experiment failed')
    raw.mkdir(parents=True,exist_ok=False);previous=json.loads((ROOT/'validation/solver_efficiency_moderate_40.json').read_text())
    defaults={}
    for mode,label in (('by','by_native'),('hispid','hi_native')):
        reference=previous['records'][mode]['after'][0]
        if digest(ROOT/reference['state'])!=reference['state_sha256']:raise RuntimeError('historical reference fingerprint changed')
        checks=[]
        for candidate in evidence['records'][label]:
            if digest(ROOT/candidate['state'])!=candidate['state_sha256']:raise RuntimeError('candidate fingerprint changed')
            checks.append(compare_states(ROOT/reference['state'],ROOT/candidate['state']))
        defaults[mode]=dict(passed=True,arrays=checks,reference_state_sha256=reference['state_sha256'])
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1');small={}
    historical=ROOT/'validation/raw/by_upstream_controls'
    for case in ('moderate_small','target_mass','bl','spin95'):
        old=json.loads((historical/(case+'_upstream.json')).read_text());oldstate=historical/(case+'_upstream.npz')
        if not old['converged'] or digest(oldstate)!=old['state_sha256']:raise RuntimeError('upstream fixture fingerprint/convergence failed')
        params=raw/(case+'_input.json');params.write_text(json.dumps(dict(config=dict(n=[12,18,8]),by_real_parameters=old['real_parameters'],by_integer_parameters=old['integer_parameters']),indent=2)+'\n')
        records={}
        for version in ('default','modal'):
            result=raw/(case+'_'+version+'.json');state=raw/(case+'_'+version+'.npz');log=raw/(case+'_'+version+'.log')
            cmd=[sys.executable,str(ROOT/'validation/compare_by_upstream.py'),'--library',args.by_library,'--parameters',str(params),'--output',str(result),'--state',str(state),'--by-preconditioner',str(int(version=='modal'))]
            if case=='target_mass':cmd+=['--target-mass']
            with log.open('w') as stream:subprocess.run(cmd,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=120)
            r=json.loads(result.read_text());r.update(command=cmd,state=str(state.relative_to(ROOT)),log=str(log.relative_to(ROOT)),log_sha256=digest(log));records[version]=r
        exact=compare_states(oldstate,ROOT/records['default']['state'])
        a=np.load(oldstate);b=np.load(ROOT/records['modal']['state']);arrays={k:dict(max_difference=float(np.max(np.abs(a[k]-b[k]))),max_scaled_difference=float(np.max(np.abs(a[k]-b[k])/(1+np.abs(a[k])))),finite=bool(np.isfinite(b[k]).all())) for k in a.files}
        charges={k:float(np.max(np.abs(np.asarray(old[k])-records['modal'][k])/(1+np.abs(old[k])))) for k in ('final_bare_masses','internal_end_adm_masses')}
        core=('v_d0','u_d0','cf_v_d0','spectral_v_samples')
        recomputed=float(np.max(np.abs(b['recomputed_F'])))
        target_ok=case!='target_mass' or max(records['modal']['target_mass_errors'])<=1e-10
        passed=all(r['converged'] for r in records.values()) and all(arrays[k]['finite'] and arrays[k]['max_scaled_difference']<=1e-10 for k in core) and all(v<=1e-10 for v in charges.values()) and target_ok and recomputed<=1e-12
        small[case]=dict(passed=bool(passed),records=records,original_default_bitwise=exact,modal_arrays=arrays,modal_charge_scaled_differences=charges,max_recomputed_F=recomputed,target_mass_gate=target_ok,upstream_state_sha256=old['state_sha256'])
        print(case,'passed',passed,'modal Newton/Krylov',records['modal']['work_statistics']['newton_iterations'],records['modal']['work_statistics']['krylov_iterations'],flush=True)
    result=dict(passed=all(r['passed'] for r in small.values()),default_full40_bitwise=defaults,small_original_controls=small,
                evidence=args.evidence,evidence_sha256=digest(args.evidence),by_library_sha256=digest(args.by_library),upstream_commit='ec563aeb672235b9443c330f9cde65f7246e8ea4',
                full_tests_log_sha256=digest(ROOT/'build-by-modal-accepted/full-tests.log'),timing_build=json.loads((ROOT/'build-by-modal-accepted/timing/build.json').read_text()),
                script_sha256={str(f.relative_to(ROOT)):digest(f) for f in (Path(__file__).resolve(),ROOT/'validation/compare_by_upstream.py')},
                note='Modal changes Krylov path, so numerical equivalence is required; inherited default remains bit-identical. Small cases test native original stopping and target-mass convergence. No high-spin binary physical acceptance inherited.')
    out.write_text(json.dumps(result,indent=2)+'\n')
    if not result['passed']:raise SystemExit(1)
if __name__=='__main__':main()
