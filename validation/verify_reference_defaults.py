"""Verify archived full states, work counts and traces after a CPU refactor."""
import argparse,json,os,subprocess,sys
from pathlib import Path
import numpy as np
from benchmark_bowen_york import ROOT,digest
from benchmark_krylov_matrix import arrays

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hispid-library',required=True);p.add_argument('--by-library',required=True);p.add_argument('--output',required=True)
    args=p.parse_args();out=Path(args.output);raw=ROOT/'validation/raw'/out.stem
    if out.exists():raise FileExistsError(out)
    raw.mkdir(parents=True,exist_ok=False)
    old=json.loads((ROOT/'validation/common_stopping_moderate_40.json').read_text());results={}
    input_path=raw/'input.json';input_path.write_text(json.dumps(old['records']['hi_native'][0]['config']))
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    for mode,label in [('hispid','hi_native'),('by','by_native')]:
        reference=old['records'][label][0];result=raw/(label+'.json');state=raw/(label+'.npz');log=raw/(label+'.log')
        command=[sys.executable,str(ROOT/'validation/benchmark_bowen_york.py'),'--hispid-library',args.hispid_library,'--by-library',args.by_library,'--case','moderate','--grid','40:80:16','--tolerance','1e-12','--output',str(out),'--mode',mode,'--worker-output',str(result),'--state-output',str(state)]
        if mode=='by':command+=['--input',str(input_path),'--by-verbose']
        with log.open('w') as stream:subprocess.run(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=1200)
        candidate=json.loads(result.read_text())
        if digest(ROOT/reference['state'])!=reference['state_sha256']:raise RuntimeError('archived reference changed')
        a=np.load(ROOT/reference['state']);b=np.load(state);metrics=arrays(a,b)
        counts=('newton_iterations','krylov_iterations') if mode=='hispid' else ('newton_iterations','krylov_iterations','jvp_applications','preconditioner_applications','relaxation_sweeps')
        refstats=reference['diagnostics'] if mode=='hispid' else reference['work_statistics'];stats=candidate['diagnostics'] if mode=='hispid' else candidate['work_statistics']
        trace=lambda path:[line for line in Path(path).read_text().splitlines() if line.startswith(('Newton:','bicgstab:'))]
        trace_equal=mode=='hispid' or trace(log)==trace(ROOT/reference['log'])
        equal_counts=all(stats[key]==refstats[key] for key in counts)
        passed=candidate['converged'] and set(a.files)==set(b.files) and all(v['bitwise_equal'] and v['finite'] for v in metrics.values()) and equal_counts and trace_equal
        results[mode]=dict(passed=bool(passed),arrays=metrics,counts_equal=equal_counts,trace_equal=trace_equal,record=candidate,state=str(state.relative_to(ROOT)),state_sha256=digest(state),log=str(log.relative_to(ROOT)),log_sha256=digest(log))
        print(mode,'passed',passed,'counts',equal_counts,'trace',trace_equal,flush=True)
    passed=all(r['passed'] for r in results.values());out.write_text(json.dumps(dict(passed=passed,controls=results),indent=2)+'\n')
    if not passed:raise SystemExit(1)
if __name__=='__main__':main()
