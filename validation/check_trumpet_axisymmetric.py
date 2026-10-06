"""Small full-system equivalence control for the optional axial Krylov sector."""
import argparse,ctypes as C,hashlib,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'examples'),str(ROOT/'validation')]
from hispid import Backend,Hole
from trumpet_configs import regularize_trumpet

def run(library,output,execution,fields_output=None):
    if output.exists():raise FileExistsError(output)
    if fields_output is not None and fields_output.exists():raise FileExistsError(fields_output)
    backend=Backend(str(library.resolve()));c=backend.config()
    c.hole[0]=Hole(.6,center=(2,0,0),velocity=(-.2,0,0))
    c.hole[1]=Hole(.4,center=(-2,0,0),velocity=(.15,0,0))
    c.omega[:]=[.5,.5];regularize_trumpet(c,12,8);c.memory_limit_mib=2048;c.krylov_restart=64
    xyz=np.array([[2.3,.4,.2],[1.7,-.5,.3],[-2.2,.3,-.4],[-1.7,-.4,-.3],[.3,1.2,.7],[4,2,-3]])
    rows=[];reference=None
    controls=[('reference',False),('reference',True)]
    if execution!='reference':controls.append((execution,True))
    for mode,axial in controls:
        start=time.monotonic()
        with backend.create(c,execution=mode,geometry='host') as s:
            if axial:
                # Deliberately contaminate the starting guess in forbidden
                # scalar and transverse-vector sectors; the accepted physical
                # solution must still match the unrestricted zero-start solve.
                guess=s.unknowns().reshape(c.n[2],c.n[1],c.n[0],4)
                guess[2,:,:,0]=1e-3
                guess[1,:,:,2]=2e-3;guess[c.n[2]//2+1,:,:,3]=-2e-3
                s.set_unknowns(guess.ravel())
            diagnostics=s.solve(krylov='gmres',linear_rtol=.1,axisymmetric=axial)
            fields=s.sample_with_derivatives(xyz)
            if reference is None:reference=fields
            errors={k:float(np.max(abs(fields[k]-reference[k])/(1+abs(reference[k])))) for k in fields}
            rows.append(dict(execution=mode,axisymmetric=axial,diagnostics=diagnostics,field_scaled_errors=errors,
                seconds=time.monotonic()-start,passed=bool(diagnostics['converged'] and max(errors.values())<1e-8)))
    # A spinning binary is outside the explicitly implemented no-swirl sector.
    c.hole[0].spin[2]=.01
    with backend.create(c) as s:
        rejected=backend.lib.HiSpID_set_axisymmetric(s.context,1)!=0
    c.hole[0].spin[2]=0;c.hole[0].velocity[1]=.01
    with backend.create(c) as s:
        transverse_rejected=backend.lib.HiSpID_set_axisymmetric(s.context,1)!=0
    result=dict(kind='axisymmetric_full_system_equivalence',rows=rows,spinning_configuration_rejected=rejected,
        transverse_boost_rejected=transverse_rejected,forbidden_initial_guess_injected=True,
        passed=bool(rejected and transverse_rejected and all(r['passed'] for r in rows)),binary_acceptance=False,
        library_sha256=backend.library_sha256(),driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if fields_output is not None:
        np.savez_compressed(fields_output,**fields)
        result['field_witness']=dict(path=str(fields_output),sha256=hashlib.sha256(fields_output.read_bytes()).hexdigest(),execution=controls[-1][0],axisymmetric=True)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    return 0 if result['passed'] else 1
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--execution',choices=['reference','kokkos'],default='reference');p.add_argument('--fields-output',type=Path);a=p.parse_args();raise SystemExit(run(a.library,a.output,a.execution,a.fields_output))
