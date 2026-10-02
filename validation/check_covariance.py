"""Solved binary rotation/translation and ADM-origin checks, one thread."""
import argparse,json,time,hashlib
from pathlib import Path
import numpy as np
from hispid import Backend,Config
from configs import moderate,as_dict
from physical import extrapolate
from run_validation import ROOT,REPORT,RAW,points
from prolong import prolong

TENSORS=('gamma','Kij','conformal_metric','Atilde')

def transformed(config,Q,c):
    new=Config.from_buffer_copy(config)
    for h in new.hole:
        h.center[:]=Q@np.array(h.center)+c
        h.spin[:]=Q@np.array(h.spin);h.velocity[:]=Q@np.array(h.velocity)
    return new

def rotate_values(values,Q):
    out={k:v.copy() for k,v in values.items()}
    for k in TENSORS:out[k]=np.einsum('ik,nkl,jl->nij',Q,values[k].reshape(-1,3,3),Q).reshape(-1,9)
    out['correction'][:,1:]=values['correction'][:,1:]@Q.T
    return out

def errors(a,b):
    out={}
    for k in (*TENSORS,'psi','mean_curvature','correction'):
        diff=np.abs(a[k]-b[k]);scale=np.maximum(np.max(np.abs(b[k]),axis=1),1e-8) if b[k].ndim==2 else np.maximum(np.abs(b[k]),1e-8)
        out[k]=float(np.max(diff/scale[:,None])) if diff.ndim==2 else float(np.max(diff/scale))
    return out

def run(backend,levels,label='moderate'):
    saved=json.loads(REPORT.read_text()).get(label,{})
    if not saved.get('passed'):raise ValueError('moderate physical convergence gate has not passed')
    def config_for(n,nphi):
        c=moderate(backend,n,nphi)
        rec=next(r for r in saved['records'] if r['resolution']==[n,n,nphi])
        c.far_radius=rec['config']['far_radius'];c.tolerance=rec['config']['tolerance']
        return c
    axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis);angle=.73
    W=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    Q=np.eye(3)+np.sin(angle)*W+(1-np.cos(angle))*(W@W);offset=np.array([.7,-.2,.4])
    records=[];base_values=[];previous_values=None;previous_shape=None
    for n,nphi in levels:
        config=config_for(n,nphi);x,near,bulk=points(config);x=x[:near+bulk]
        start=time.monotonic()
        with backend.create(config) as s:
            s.set_unknowns(np.load(RAW/f'{label}_{n}_{nphi}.npz')['unknowns'])
            original=s.sample(x);base_values.append(original)
            prime=transformed(config,Q,offset)
            with backend.create(prime) as t:
                if previous_values is not None:t.set_unknowns(prolong(previous_values,previous_shape,list(prime.n)))
                diagnostic=t.solve();rotated=t.sample(x@Q.T+offset)
                previous_values=t.unknowns();previous_shape=list(prime.n)
                err=errors(rotated,rotate_values(original,Q))
                radii=[100.,200.,400.]
                qrot=extrapolate(radii,[t.charges(R,center=offset,ntheta=12,nphi=24) for R in radii])
                qglobal=extrapolate(radii,[t.charges(R,ntheta=12,nphi=24) for R in radii])
            base_record=next(r for r in saved['records'] if r['resolution']==[n,n,nphi])
            base_charge=np.array(base_record['charges_extrapolated']);expected=np.r_[base_charge[0],Q@base_charge[1:4],Q@base_charge[4:]]
            origin_expected=expected.copy();origin_expected[4:]+=np.cross(offset,expected[1:4])
            RAW.mkdir(exist_ok=True);np.savez_compressed(RAW/f'covariance_{n}_{nphi}.npz',points=x,**{'base_'+k:v for k,v in original.items()},**{'rotated_'+k:v for k,v in rotated.items()})
        rec=dict(resolution=[n,n,nphi],library_sha256=hashlib.sha256(backend.path.read_bytes()).hexdigest(),rotation=Q.tolist(),offset=offset.tolist(),config=as_dict(prime),diagnostics=diagnostic,
                 relative_errors=err,rotated_EPJ=qrot.tolist(),global_origin_EPJ=qglobal.tolist(),expected_EPJ=expected.tolist(),
                 charge_error=np.abs(qrot-expected).tolist(),origin_charge_error=np.abs(qglobal-origin_expected).tolist(),seconds=time.monotonic()-start)
        records.append(rec);print(json.dumps(rec),flush=True)
        report=json.loads(REPORT.read_text());report['covariance']=dict(records=records,passed=False);REPORT.write_text(json.dumps(report,indent=2)+'\n')
    truncation=errors(base_values[-2],base_values[-1]);last=records[-1]
    field_pass=all(last['relative_errors'][k]<=max(1e-9,5*truncation[k]) for k in truncation)
    charge_pass=max(last['charge_error']+last['origin_charge_error'])<.005
    # Translation alone leaves the canonical-frame coefficients invariant.
    n,nphi=levels[-1];cfg=config_for(n,nphi);x,near,bulk=points(cfg);x=x[:near+bulk]
    with backend.create(transformed(cfg,np.eye(3),offset)) as t:
        t.set_unknowns(np.load(RAW/f'{label}_{n}_{nphi}.npz')['unknowns']);d=t.solve()
        translation_errors=errors(t.sample(x+offset),base_values[-1])
    translation_pass=max(translation_errors.values())<1e-9 and d['status']==0
    return dict(source_case=label,records=records,finest_base_truncation=truncation,translation_errors=translation_errors,translation_diagnostics=d,
                field_pass=bool(field_pass),charge_pass=bool(charge_pass),translation_pass=bool(translation_pass),
                passed=bool(field_pass and charge_pass and translation_pass and all(r['diagnostics']['status']==0 for r in records)))

def main():
    p=argparse.ArgumentParser();p.add_argument('--library',required=True);p.add_argument('--levels',default='24:12,40:20,56:28');p.add_argument('--case',default='moderate');a=p.parse_args()
    b=Backend(a.library);r=run(b,[tuple(map(int,s.split(':'))) for s in a.levels.split(',')],a.case)
    report=json.loads(REPORT.read_text());report['covariance']=r;REPORT.write_text(json.dumps(report,indent=2)+'\n')
    return 0 if r['passed'] else 1

if __name__=='__main__':raise SystemExit(main())
