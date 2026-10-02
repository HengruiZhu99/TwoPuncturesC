"""Refine ADM quadrature on a physically checked convergence sequence.

No solves, checkpoint relabeling, or binary horizon claims are performed.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from hispid import Backend
from checkpoints import ROOT,restore_equivalent
from physical import extrapolate

def main():
    p=argparse.ArgumentParser();p.add_argument('--library',required=True)
    p.add_argument('--proof',required=True);p.add_argument('--case',required=True)
    p.add_argument('--output',required=True);p.add_argument('--radii',default='100,200,400')
    p.add_argument('--angular-bound',type=float,default=1e-7)
    a=p.parse_args();backend=Backend(a.library);proof=json.loads(Path(a.proof).read_text())
    source=json.loads((ROOT/'validation/results.json').read_text())[a.case]
    records=source['records'][-3:];radii=list(map(float,a.radii.split(',')))
    if len(records)!=3 or len(radii)<3 or not np.all(np.isfinite(radii)) or not np.all(np.diff(radii)>0):
        raise ValueError('three spectral levels and >=3 increasing radii required')
    if not 0<a.angular_bound<=1e-7:raise ValueError('positive angular bound <=1e-7 required')
    result=dict(case=a.case,library_sha256=backend.library_sha256(),compatibility_proof=a.proof,
        radii=radii,angular_absolute_bound=a.angular_bound,spectral_charge_bound=.005,
        records=[],preliminary_physical_and_charge_gate=False,stronger_gate=False,
        covariance_verified=False,horizon_enclosure_verified=False,binary_validation_complete=False)
    fits=[]
    for rec in records:
        cfg,values=restore_equivalent(backend,rec,proof);nb=rec['resolution'][1]
        row=dict(resolution=rec['resolution'],checkpoint_source_sha256=rec['library_sha256'],quadratures=[])
        start=time.monotonic()
        with backend.create_sampler(cfg) as solution:
            solution.set_unknowns(values)
            for nt,np_ in [(2*nb,64),(3*nb,64),(3*nb,128)]:
                before=time.monotonic();q=np.array([solution.charges(r,ntheta=nt,nphi=np_) for r in radii])
                row['quadratures'].append(dict(ntheta=nt,nphi=np_,EPJ=q.tolist(),
                    extrapolated_EPJ=extrapolate(radii,q).tolist(),seconds=time.monotonic()-before))
                Path(a.output).write_text(json.dumps(result|{'in_progress':row},indent=2)+'\n')
                print(rec['resolution'],nt,np_,row['quadratures'][-1],flush=True)
        q=[np.array(v['EPJ']) for v in row['quadratures']]
        f=[np.array(v['extrapolated_EPJ']) for v in row['quadratures']]
        row['polar_charge_change_linf']=float(np.max(abs(q[1]-q[0])))
        row['azimuthal_charge_change_linf']=float(np.max(abs(q[2]-q[1])))
        row['extrapolated_quadrature_change_linf']=float(max(np.max(abs(f[1]-f[0])),np.max(abs(f[2]-f[1]))))
        row['angular_verified']=bool(max(row['polar_charge_change_linf'],row['azimuthal_charge_change_linf'],row['extrapolated_quadrature_change_linf'])<a.angular_bound)
        row['seconds']=time.monotonic()-start;result['records'].append(row);fits.append(f[-1])
        Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    result['angular_verified']=all(r['angular_verified'] for r in result['records'])
    result['finest_spectral_charge_change_linf']=float(np.max(abs(fits[-1]-fits[-2])))
    result['spectral_charge_stable']=result['finest_spectral_charge_change_linf']<.005
    result['physical_sequence_improves']=all(records[i+1][k][q]<records[i][k][q] for i in range(2) for k in ('near','bulk') for q in ('H_rms','M_rms'))
    best=records[-1]
    result['preliminary_physical_and_charge_gate']=bool(best['passed_local'] and result['physical_sequence_improves'] and result['angular_verified'] and result['spectral_charge_stable'])
    result['stronger_gate']=bool(best['passed_strict'] and result['preliminary_physical_and_charge_gate'])
    result['note']='Read-only exact-SHA field/operator migration verifies integration on the original checked geometries. Coarse charge summaries and original failed records remain unchanged. Covariance, stronger accuracy and binary horizons/enclosure remain separate requirements.'
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
    if not result['preliminary_physical_and_charge_gate']:raise SystemExit(1)

if __name__=='__main__':main()
