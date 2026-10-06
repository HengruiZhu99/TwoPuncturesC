"""Measure the no-swirl axisymmetric modal sector of coaxial, nonspinning data.

No coefficients are changed and no solution acceptance is inferred. Scalar/u
and b_x have only m=0; b_y and b_z share the cos(phi)/sin(phi) coefficient.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def inspect(run):
    meta=json.loads((run/'result.json').read_text());c=meta['config']
    for h in c['hole']:
        if any(h['spin']) or any(h['velocity'][1:]) or any(h['center'][1:]):
            raise ValueError('coaxial nonspinning configuration required')
    na,nb,nphi=c['n']
    with np.load(run/'solve.npz') as saved:
        v=saved['unknowns'].reshape(nphi,nb,na,4)
    projected=np.zeros_like(v)
    projected[0,:,:,0:2]=v[0,:,:,0:2]
    radial=.5*(v[1,:,:,2]+v[nphi//2+1,:,:,3])
    projected[1,:,:,2]=radial;projected[nphi//2+1,:,:,3]=radial
    amplitude=np.max(abs(v),axis=(0,1,2));error=np.max(abs(v-projected),axis=(0,1,2))
    return dict(run=str(run),resolution=c['n'],maximum_nodal_by_component=amplitude.tolist(),
        forbidden_sector_max_by_component=error.tolist(),
        forbidden_sector_relative_by_component=(error/(1+amplitude)).tolist(),
        maximum_nodal_by_mode_component=np.max(abs(v),axis=(1,2)).tolist(),
        result_sha256=hashlib.sha256((run/'result.json').read_bytes()).hexdigest(),
        solve_sha256=hashlib.sha256((run/'solve.npz').read_bytes()).hexdigest())

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    result=dict(rows=[inspect(run) for run in a.runs],binary_acceptance=False,
        note='Auxiliary modal P coefficients only, not physical-field errors or a proof that symmetry restriction fixes the solve.',
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    for row in result['rows']:print(row['run'],row['forbidden_sector_relative_by_component'])
