"""Compare FD/analytic native ADM integrals in image-verified processes."""
import argparse,hashlib,json,os,subprocess,sys,time
from pathlib import Path
import numpy as np
from hispid import Backend
from checkpoints import ROOT,select_record,restore_payload

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--old-library');p.add_argument('--new-library')
    p.add_argument('--equivalence-proof');p.add_argument('--output',required=True)
    p.add_argument('--case',default='moderate_default_polar2_r128');p.add_argument('--radii',default='40,200,1000')
    p.add_argument('--quadrature',default='16:32');p.add_argument('--worker-library');p.add_argument('--source-sha')
    p.add_argument('--momentum-atol',type=float,default=0,help='predeclared P/J rounding bound; default requires bitwise equality')
    a=p.parse_args();record=select_record(a.case);nt,np_=map(int,a.quadrature.split(':'))
    if not np.isfinite(a.momentum_atol) or a.momentum_atol<0:raise ValueError('finite nonnegative P/J bound required')
    if a.worker_library:
        if record['library_sha256']!=a.source_sha:raise ValueError('wrong source checkpoint')
        backend=Backend(a.worker_library)
        if backend.parameterization()!=record['unknown_parameterization_id'] or backend.parameterization_maps()!=record['collocation_maps']:
            raise ValueError('different continuous basis/maps')
        cfg,values=restore_payload(record,backend.config())
        start=time.monotonic()
        with backend.create_sampler(cfg) as solution:
            solution.set_unknowns(values)
            charges=[solution.charges(float(r),ntheta=nt,nphi=np_).tolist() for r in a.radii.split(',')]
        Path(a.output).write_text(json.dumps(dict(library_sha256=backend.library_sha256(),charges=charges,seconds=time.monotonic()-start),indent=2)+'\n')
        return
    proof=json.loads(Path(a.equivalence_proof).read_text())
    hashes=[sha(path) for path in (a.old_library,a.new_library)]
    if hashes!=[proof['old_library_sha256'],proof['new_library_sha256']] or not proof['passed_off_axis_and_equations'] or not proof['equations_bitwise_identical']:
        raise ValueError('matching bitwise field/operator proof required before method comparison')
    temporary=ROOT/'validation/raw'/Path(a.output).stem;temporary.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    reports=[]
    for side,library in zip(('old','new'),(a.old_library,a.new_library)):
        path=temporary/(side+'.json')
        command=[sys.executable,str(Path(__file__).resolve()),'--worker-library',str(Path(library).resolve()),'--source-sha',hashes[0],
                 '--case',a.case,'--radii',a.radii,'--quadrature',a.quadrature,'--output',str(path)]
        subprocess.run(command,cwd=ROOT,env=env,check=True,timeout=1200)
        reports.append(json.loads(path.read_text()))
    if [r['library_sha256'] for r in reports]!=hashes:raise ValueError('worker image mismatch')
    difference=np.abs(np.array(reports[0]['charges'])-reports[1]['charges'])
    result=dict(case=a.case,resolution=record['resolution'],radii=list(map(float,a.radii.split(','))),quadrature=[nt,np_],
        image_verified_separate_processes=True,records=reports,absolute_charge_differences=difference.tolist(),
        energy_absolute_tolerance=1e-8,momentum_and_angular_momentum_require_bitwise=a.momentum_atol==0,momentum_and_angular_momentum_absolute_tolerance=a.momentum_atol,
        speedup=reports[0]['seconds']/reports[1]['seconds'],
        passed=bool(np.max(difference[:,0])<1e-8 and np.max(difference[:,1:])<=a.momentum_atol),
        binary_acceptance=False,note='Charge integration implementation comparison at identical declared sphere nodes; loaded images and source payload are verified independently. Failed physical binary flags are unchanged.')
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if not result['passed']:raise SystemExit(1)

if __name__=='__main__':main()
