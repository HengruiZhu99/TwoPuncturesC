"""Serial cross-backend BY original fixtures, including target-mass iteration."""
import argparse,json,os,subprocess,sys
from pathlib import Path
import numpy as np
from benchmark_bowen_york import ROOT,digest


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',required=True);p.add_argument('--output',required=True);p.add_argument('--execution',choices=('reference','kokkos'),default='reference');args=p.parse_args()
    out=Path(args.output);raw=ROOT/'validation/raw'/out.stem
    if out.exists():raise FileExistsError(out)
    raw.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    records={}
    for case in ('moderate_small','target_mass','bl','spin95'):
        historical=ROOT/'validation/raw/by_upstream_controls';old=json.loads((historical/(case+'_upstream.json')).read_text())
        oldstate=historical/(case+'_upstream.npz')
        if digest(oldstate)!=old['state_sha256']:raise RuntimeError('original fixture hash changed')
        a=np.load(oldstate);params=raw/(case+'_input.json')
        params.write_text(json.dumps(dict(config=dict(n=[12,18,8]),by_real_parameters=old['real_parameters'],by_integer_parameters=old['integer_parameters']),indent=2)+'\n')
        series={}
        for krylov in ('bicgstab','gmres'):
            label=case+'_'+krylov;result=raw/(label+'.json');state=raw/(label+'.npz');log=raw/(label+'.log')
            cmd=[sys.executable,str(ROOT/'validation/compare_by_upstream.py'),'--library',args.library,'--parameters',str(params),
                 '--output',str(result),'--state',str(state),'--by-preconditioner','1','--krylov',krylov,'--execution',args.execution]
            if case=='target_mass':cmd+=['--target-mass']
            with log.open('w') as stream:subprocess.run(cmd,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=180,check=True)
            r=json.loads(result.read_text());b=np.load(state)
            metrics={k:dict(finite=bool(np.isfinite(b[k]).all()),scaled_difference=float(np.max(np.abs(a[k]-b[k])/(1+np.abs(a[k]))))) for k in a.files}
            charges={k:float(np.max(np.abs(np.asarray(old[k])-r[k])/(1+np.abs(old[k])))) for k in ('final_bare_masses','internal_end_adm_masses')}
            passed=r['converged'] and all(v['finite'] for v in metrics.values()) and all(metrics[k]['scaled_difference']<=1e-10 for k in ('v_d0','u_d0','cf_v_d0','spectral_v_samples')) and all(v<=1e-10 for v in charges.values()) and float(np.max(np.abs(b['recomputed_F'])))<=1e-12
            if case=='target_mass':passed=passed and max(r['target_mass_errors'])<=1e-10
            r.update(passed=bool(passed),arrays=metrics,charges_scaled_difference=charges,state=str(state.relative_to(ROOT)),log=str(log.relative_to(ROOT)),log_sha256=digest(log),command=cmd)
            series[krylov]=r;print(label,'passed',passed,flush=True)
        records[case]=series
    result=dict(passed=all(r['passed'] for series in records.values() for r in series.values()),records=records,execution=args.execution,library_sha256=digest(args.library),script_sha256=digest(__file__),original_commit='ec563aeb672235b9443c330f9cde65f7246e8ea4',note='Modal M with both methods; opt-in Kokkos explicitly uses true RHS-relative inner tolerance. Native outer tolerance retained; not physical binary acceptance.')
    out.write_text(json.dumps(result,indent=2)+'\n')
    if not result['passed']:raise SystemExit(1)
if __name__=='__main__':main()
