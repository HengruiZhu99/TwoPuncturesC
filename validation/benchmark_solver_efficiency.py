"""Fresh-process, serial before/after timing and exact state controls.

One numerical worker at a time. Same backend configurations/stopping criteria;
BY/HiSpID solutions are never compared for physical equivalence.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import numpy as np
from benchmark_bowen_york import ROOT, digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('hispid-before','hispid-after','by-before','by-after'):
        p.add_argument('--'+name,required=True)
    p.add_argument('--grid',default='40:80:16')
    p.add_argument('--case',default='moderate',choices=('moderate','spin95','boost885'))
    p.add_argument('--tolerance',type=float,default=1e-12)
    p.add_argument('--repeats',type=int,default=2)
    p.add_argument('--output',required=True)
    args=p.parse_args();out=Path(args.output);raw=ROOT/'validation/raw'/out.stem
    if out.exists():raise FileExistsError(out)
    raw.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    records={};input_path=raw/'input.json'
    for mode in ('hispid','by'):
        records[mode]={'before':[],'after':[]}
        for repeat in range(args.repeats):
            for version in (('before','after') if repeat%2==0 else ('after','before')):
                label=f'{mode}_{version}_{repeat}';state=raw/(label+'.npz');result=raw/(label+'.json');log=raw/(label+'.log')
                cmd=[sys.executable,str(ROOT/'validation/benchmark_bowen_york.py'),
                     '--hispid-library',getattr(args,'hispid_'+version),
                     '--by-library',getattr(args,'by_'+version),'--case',args.case,
                     '--grid',args.grid,'--tolerance',str(args.tolerance),'--output',str(out),
                     '--mode',mode,'--worker-output',str(result),'--state-output',str(state)]
                if mode=='by':cmd+=['--input',str(input_path),'--by-verbose']
                print('Starting',label,flush=True)
                with log.open('w') as stream:subprocess.run(cmd,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=1200)
                r=json.loads(result.read_text());r.update(state=str(state.relative_to(ROOT)),state_sha256=digest(state),log=str(log.relative_to(ROOT)),log_sha256=digest(log),command=cmd)
                records[mode][version].append(r)
                if mode=='hispid' and not input_path.exists():input_path.write_text(json.dumps(r['config'],indent=2)+'\n')
                print(label,round(r['solve_seconds'],6),r['max_rss_bytes'],r['converged'],flush=True)
    invariance={};performance={}
    for mode,versions in records.items():
        reference=versions['before'][0];state0=np.load(ROOT/reference['state']);checks=[]
        for version,series in versions.items():
            for r in series:
                state=np.load(ROOT/r['state']);same_keys=set(state.files)==set(state0.files) and all(state[k].shape==state0[k].shape and state[k].dtype==state0[k].dtype for k in state0.files)
                arrays={name:dict(bitwise_equal=state[name].tobytes()==state0[name].tobytes(),
                                 max_absolute_difference=float(np.max(np.abs(state[name]-state0[name]))),
                                 finite=bool(np.isfinite(state[name]).all()),elements=state[name].size)
                        for name in state0.files} if same_keys else {}
                trace=lambda record:[line for line in (ROOT/record['log']).read_text().splitlines() if line.startswith(('Newton:','bicgstab:'))]
                parameters_equal=r['config']==reference['config'] and all(r.get(k)==reference.get(k) for k in ('by_real_parameters','by_integer_parameters'))
                diagnostics_equal={k:v for k,v in r['diagnostics'].items() if k!='seconds'}=={k:v for k,v in reference['diagnostics'].items() if k!='seconds'}
                checks.append(dict(parameters_equal=parameters_equal,diagnostics_equal=diagnostics_equal,version=version,converged=r['converged'],same_arrays=same_keys,arrays=arrays,
                                   by_iteration_trace_bitwise_equal=bool(trace(r)) and trace(r)==trace(reference) if mode=='by' else None,
                                   hispid_iteration_counts_equal=all(r['diagnostics'][k]==reference['diagnostics'][k] for k in ('newton_iterations','krylov_iterations')) if mode=='hispid' else None))
        passed=all(c['converged'] and c['parameters_equal'] and c['diagnostics_equal'] and c['same_arrays'] and all(a['bitwise_equal'] and a['finite'] for a in c['arrays'].values()) and c['by_iteration_trace_bitwise_equal'] is not False and c['hispid_iteration_counts_equal'] is not False for c in checks)
        invariance[mode]=dict(passed=passed,checks=checks)
        medians={version:{key:float(np.median([r[key] for r in series])) for key in ('solve_seconds','ready_to_sample_seconds','max_rss_bytes')} for version,series in versions.items()}
        performance[mode]=dict(medians=medians,
                              solve_speedup=medians['before']['solve_seconds']/medians['after']['solve_seconds'] if passed else None,
                              ready_speedup=medians['before']['ready_to_sample_seconds']/medians['after']['ready_to_sample_seconds'] if passed else None,
                              peak_rss_reduction_fraction=1-medians['after']['max_rss_bytes']/medians['before']['max_rss_bytes'] if passed else None)
    result=dict(case=args.case,grid=list(map(int,args.grid.split(':'))),tolerance=args.tolerance,repeats=args.repeats,
                records=records,invariance=invariance,performance=performance,
                cpu_threads=1,isolated_sequential_workers=True,platform=platform.platform(),python=sys.version,
                note='Native weights differ: BY sin^3(alpha)sin^3(beta), HiSpID sin^6(alpha)sin^6(beta). Both criteria unchanged before/after. State capture occurs after timing/RSS. Full-precision BY stdout trace is enabled identically in both builds. These tests do not inherit binary physical acceptance.')
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(performance,indent=2),flush=True)
    if not all(r['passed'] for r in invariance.values()):raise SystemExit(1)
if __name__=='__main__':main()
