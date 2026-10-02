"""Values-only reconstruction diagnostic for failed Hi auxiliary-field gates.

Read retained states; no native library or new solve. Preserve strict failure.
"""
import argparse,json
from pathlib import Path
import numpy as np
from benchmark_bowen_york import ROOT,digest


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',required=True);p.add_argument('--output',required=True);args=p.parse_args()
    out=Path(args.output)
    if out.exists():raise FileExistsError(out)
    result=json.loads(Path(args.evidence).read_text());a_record=result['records']['hi_gmres'][0];b_record=result['records']['hi_bicgstab'][0]
    for record in (a_record,b_record):
        if digest(ROOT/record['state'])!=record['state_sha256']:raise RuntimeError('retained state changed')
    na,nb,nphi=a_record['config']['n'];a=np.load(ROOT/a_record['state'])['unknowns'].reshape(nphi,nb,na,4);b=np.load(ROOT/b_record['state'])['unknowns'].reshape(a.shape)
    # These runs use the explicit default .2/2 maps in HiSpID_axis.hpp.
    s=.5*(1-np.cos(np.pi*(np.arange(na)+.5)/na));t=.2*s/(1-.8*s)
    eta=np.tanh(-2*np.cos(np.pi*(np.arange(nb)+.5)/nb))/np.tanh(2)
    q=np.sqrt(t[None,:]*(1-eta[:,None]**2));factor=np.empty((nphi,nb,na))
    phi=2*np.pi*np.arange(nphi)/nphi;inverse=np.empty((nphi,nphi));half=nphi//2
    for k in range(nphi):
        mode=k if k<=half else k-half;r=mode if mode<=4 else 3 if mode%2 else 4
        factor[k]=-2*(1-t[None,:])*q**r
        inverse[:,k]=np.sqrt((1 if mode in (0,half) else 2)/nphi)*(np.cos(mode*phi) if k<=half else np.sin(mode*phi))
    modal_delta=(b-a)*factor[:,:,:,None]
    reconstructed=np.einsum('lk,kjiv->ljiv',inverse,modal_delta,optimize=True)
    scaled=np.abs(b-a)/(1+np.abs(a));index=tuple(map(int,np.unravel_index(np.argmax(scaled),a.shape)))
    r=dict(evidence=args.evidence,evidence_sha256=digest(args.evidence),script_sha256=digest(__file__),
      unknown_gate_passed=False,coefficient_gate_passed=False,max_deltaP_node=dict(mode_index_B_A_component=index,reference=float(a[index]),candidate=float(b[index]),scaled_difference=float(scaled[index]),physical_modal_factor=float(factor[index[:3]]),physical_modal_delta=float(modal_delta[index])),
      modal_max_scaled_differences=[float(np.max(scaled[k])) for k in range(nphi)],
      maximum_physical_correction_value_difference_over_all_collocation_nodes=np.max(np.abs(reconstructed),axis=(0,1,2)).tolist(),
      limits='Values-only algebraic reconstruction with fixed maps/basis. No gradient/Hessian, off-grid or exact-nullspace certification. Does not waive failed1e-10 auxiliary/coefficient gates.')
    out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
