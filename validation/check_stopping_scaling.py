"""Fresh-process positive row-scaling control with unchanged physical fields."""
import argparse,json,os,subprocess,sys
from pathlib import Path
import numpy as np
from benchmark_bowen_york import ROOT,configuration,digest,SAMPLE_POINTS
from hispid import Backend


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('hi6','hi3','archived6','output'):p.add_argument('--'+k,required=True)
    p.add_argument('--worker');p.add_argument('--state');args=p.parse_args()
    if args.worker:
        b=Backend(args.worker);c=configuration(b,'moderate',(12,18,8),1e-12)
        with b.create(c) as data:
            rng=np.random.default_rng(14108607);v=rng.normal(0,1e-6,data.size);direction=rng.normal(0,1e-4,data.size)
            data.set_unknowns(v)
            np.savez(args.state,residual=data.residual(v),jvp=data.jvp(v,direction),**data.sample_with_derivatives(SAMPLE_POINTS))
        return
    out=Path(args.output);raw=ROOT/'validation/raw'/out.stem
    if out.exists():raise FileExistsError(out)
    raw.mkdir(parents=True,exist_ok=False);files={};commands={}
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    for label in ('archived6','hi6','hi3'):
        file=raw/(label+'.npz');cmd=[sys.executable,str(Path(__file__).resolve()),'--hi6',args.hi6,'--hi3',args.hi3,'--archived6',args.archived6,'--output',str(out),'--worker',getattr(args,label),'--state',str(file)]
        with (raw/(label+'.log')).open('w') as stream:subprocess.run(cmd,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=120)
        files[label]=dict(path=str(file.relative_to(ROOT)),sha256=digest(file),library_sha256=digest(getattr(args,label)));commands[label]=cmd
    old=np.load(ROOT/files['archived6']['path']);six=np.load(ROOT/files['hi6']['path']);three=np.load(ROOT/files['hi3']['path'])
    unchanged={k:dict(bitwise_equal=old[k].tobytes()==six[k].tobytes(),max_difference=float(np.max(np.abs(old[k]-six[k])))) for k in old.files}
    fields={k:dict(bitwise_equal=six[k].tobytes()==three[k].tobytes(),max_difference=float(np.max(np.abs(six[k]-three[k])))) for k in six.files if k not in ('residual','jvp')}
    sn=np.broadcast_to(np.sin(np.pi*(np.arange(12)+.5)/12)[None,None,:]*np.sin(np.pi*(np.arange(18)+.5)/18)[None,:,None],(8,18,12)).ravel()
    scaling={}
    for k in ('residual','jvp'):
        expected=six[k].reshape(-1,4)/sn[:,None]**3;actual=three[k].reshape(-1,4)
        delta=float(np.max(np.abs(expected-actual)/(1+np.abs(expected))))
        scaling[k]=dict(max_scaled_difference=delta,bound=2e-13,passed=delta<=2e-13)
    passed=all(r['bitwise_equal'] for r in unchanged.values()) and all(r['bitwise_equal'] for r in fields.values()) and all(r['passed'] for r in scaling.values())
    result=dict(passed=bool(passed),files=files,commands=commands,default6_vs_archived_bitwise=unchanged,cubic_vs_six_physical_fields=fields,positive_scaling=scaling,grid=[12,18,8],fixture_seed=14108607)
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(scaling,indent=2),flush=True)
    if not passed:raise SystemExit(1)
if __name__=='__main__':main()
