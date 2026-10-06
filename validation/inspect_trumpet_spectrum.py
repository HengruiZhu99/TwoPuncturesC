"""Diagnose retained meridional polynomial tails without solving or acceptance."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np

def transform(values, axis):
    n=values.shape[axis]
    doubled=np.concatenate((values,np.flip(values,axis=axis)),axis=axis)
    f=np.fft.rfft(doubled,axis=axis)
    shape=[1]*values.ndim;shape[axis]=n
    phase=(np.exp(-.5j*np.pi*np.arange(n)/n)*(-1.)**np.arange(n)).reshape(shape)
    return (np.take(f,np.arange(n),axis=axis)*phase).real/n

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    rows=[]
    for run in a.runs:
        meta=json.loads((run/'result.json').read_text());na,nb,nphi=meta['config']['n']
        with np.load(run/'solve.npz') as z:values=z['unknowns'].reshape(nphi,nb,na,4)
        if not np.isfinite(values).all():raise ValueError('nonfinite nodal coefficients')
        coefficients=transform(transform(values,2),1)
        if not np.isfinite(coefficients).all():raise ValueError('nonfinite transformed coefficients')
        # Check the FFT transform convention against direct native Chebyshev
        # coefficient sums on a small selection of actual input lines.
        witness=0.
        for axis,n in ((2,na),(1,nb)):
            moved=np.moveaxis(values,axis,-1).reshape(-1,n)
            lines=moved[[0,len(moved)//2,-1]]
            modes=np.arange(n);nodes=np.arange(n)+.5
            matrix=(2./n)*(-1.)**modes[:,None]*np.cos(np.pi*modes[:,None]*nodes/n)
            direct=np.einsum('ij,kj->ik',lines,matrix,optimize=False)
            if not np.isfinite(direct).all():raise ValueError('nonfinite direct coefficient witness')
            witness=max(witness,float(np.max(abs(transform(lines,1)-direct))/(1+np.max(abs(direct)))))
        if witness>1e-11:raise ValueError('FFT/direct coefficient convention mismatch')
        axes={}
        for name,axis,n in (('radial',2,na),('polar',1,nb)):
            envelope=np.max(abs(coefficients),axis=tuple(j for j in range(3) if j!=axis))
            tail=np.take(coefficients,np.arange(7*n//8,n),axis=axis)
            axes[name]=dict(final_eighth_max_by_mode_component=np.max(abs(tail),axis=(1,2)).tolist(),envelope_by_degree_component=envelope.tolist(),
                tail_quarter_max_by_component=np.max(envelope[3*n//4:],axis=0).tolist(),
                final_eighth_max_by_component=np.max(envelope[7*n//8:],axis=0).tolist())
        rows.append(dict(resolution=[na,nb,nphi],axes=axes,fft_direct_scaled_error=witness,
            maximum_nodal_by_component=np.max(abs(values),axis=(0,1,2)).tolist(),
            solve_sha256=hashlib.sha256((run/'solve.npz').read_bytes()).hexdigest(),
            checkpoint=meta['checkpoint']))
    a.output.write_text(json.dumps(dict(rows=rows,binary_acceptance=False,
        purpose='Polynomial-tail diagnosis only; modal P tails are not physical constraint norms',
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
    for row in rows:print(row['resolution'],{k:v['final_eighth_max_by_component'] for k,v in row['axes'].items()},flush=True)

if __name__=='__main__':main()
