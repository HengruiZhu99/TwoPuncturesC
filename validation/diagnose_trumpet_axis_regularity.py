"""Inspect the forbidden quartic m6 axis term of retained C2 scalar polynomials.

Away from punctures, a smooth Cartesian scalar has Fourier mode m6=O(rho^6).
The current capped basis allows rho^4 P6. Its quartic coefficient is
-(1-t)^5 P6/(8 b^4), t=a^2, so P6 must vanish on each axis boundary.
This reports a representation diagnostic, not a physical acceptance test.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np

def weights(n, target):
    theta=np.pi*(np.arange(n)+.5)/n
    nodes=-np.cos(theta)
    hit=np.flatnonzero(abs(target-nodes)<2e-15)
    if len(hit):
        out=np.zeros(n);out[hit[0]]=1;return out
    w=(-1.)**np.arange(n)*np.sin(theta)/(target-nodes)
    return w/np.sum(w)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    rows=[]
    for run in a.runs:
        record=json.loads((run/'result.json').read_text());na,nb,nphi=record['config']['n']
        if nphi!=16:raise ValueError('this probe targets the retained moderate nphi16 sequence')
        with np.load(run/'solve.npz') as z:values=z['unknowns'].reshape(nphi,nb,na,4)[[6,14],:,:,0].copy()
        # This retained sequence uses the original .2/2 map. Verify metadata.
        with (run/'diagnostic.checkpoint').open() as stream:
            signature=stream.readline().strip();parameterization=stream.readline().strip()
        if signature!='HISPID_CHECKPOINT 2' or parameterization!='parameterization modal_P_C2prolate_mapped_v2':raise ValueError('expected v2 checkpoint with original map')
        if hashlib.sha256((run/'solve.npz').read_bytes()).hexdigest()!=record['solve_artifact_sha256']:raise ValueError('retained vector hash mismatch')
        b=np.linalg.norm(np.subtract(record['config']['hole'][0]['center'],record['config']['hole'][1]['center']))/2
        samples=[]
        for t in (.04,.25,.64):
            sigma=t/(.2+.8*t);wa=weights(na,2*sigma-1)
            for eta in (-1.,1.):
                wb=weights(nb,eta);P=np.einsum('kji,j,i->k',values,wb,wa,optimize=False)
                coeff=-(1-t)**5*P/(8*b**4)
                samples.append(dict(axis='exterior',t=t,eta=eta,x=eta*b*(1+t)/(1-t),P6=P.tolist(),quartic_coefficients=coeff.tolist()))
        for eta in (-.75,0.,.75):
            zeta=np.arctanh(eta*np.tanh(2.))/2
            P=np.einsum('kji,j,i->k',values,weights(nb,zeta),weights(na,-1.),optimize=False)
            samples.append(dict(axis='between_holes',t=0.,eta=eta,x=b*eta,P6=P.tolist(),quartic_coefficients=(-P/(8*b**4)).tolist()))
        if not np.isfinite([s['quartic_coefficients'] for s in samples]).all():raise ValueError('nonfinite axis trace')
        rows.append(dict(resolution=[na,nb,nphi],checkpoint=record['checkpoint'],solve_sha256=hashlib.sha256((run/'solve.npz').read_bytes()).hexdigest(),samples=samples))
    a.output.write_text(json.dumps(dict(purpose=__doc__,rows=rows,binary_acceptance=False,driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
    for row in rows:print(row['resolution'],max(np.linalg.norm(s['quartic_coefficients']) for s in row['samples']),flush=True)
if __name__=='__main__':main()
