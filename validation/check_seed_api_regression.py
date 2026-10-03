"""Bitwise replay of default host seed/operator/sampler API values.

Capture once from the pre-change image, then replay in a fresh process using
the replacement image. This is an API/math regression, not a solved binary.
Both the fixture and native image/dependency hashes are retained separately.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from hispid import Backend,Hole


def snapshot(backend):
    result={}
    xyz=np.array([[.125,0,0],[0,.25,0],[0,0,.5],[.1,.2,.3],[1,-2,3],[100,200,-300]])
    holes=[Hole(1),Hole(1,spin=(.57,.76,0)),Hole(1,velocity=(.531,.708,0)),
           Hole(1,spin=(.2,.3,.4),velocity=(.3,-.2,.1)),Hole(1,spin=(0,0,.99)),
           Hole(1,velocity=(np.sqrt(.99),0,0))]
    for kind,hole in enumerate(holes):
        for choice in range(2):
            for field,value in backend.seed(hole,xyz,choice).items():
                result[f'seed_{kind}_{choice}_{field}']=value
    config=backend.config();config.n[:]=[6,6,4]
    config.hole[0]=Hole(.6,(3,0,0),(.072,.054,.108),(.03,.06,.01))
    config.hole[1]=Hole(.4,(-3,0,0),(-.032,.048,.016),(-.02,-.07,.025))
    with backend.create_sampler(config) as solution:
        solution.set_unknowns(1e-7*np.cos(.23*np.arange(solution.size)))
        for field,value in solution.sample_with_derivatives(xyz+(.9,.6,-.4)).items():
            result['sampler_'+field]=value
    result['operators']=np.array([backend.operators(config,x,1e-6*np.sin(.31*np.arange(40))) for x in xyz+(.9,.6,-.4)])
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library',required=True);p.add_argument('--baseline',required=True)
    p.add_argument('--capture',action='store_true');p.add_argument('--receipt',required=True)
    a=p.parse_args();backend=Backend(a.library);values=snapshot(backend)
    images={str(backend.path):backend.library_sha256(),**backend.dependency_images}
    baseline=Path(a.baseline);receipt=Path(a.receipt)
    if a.capture:
        baseline.parent.mkdir(parents=True,exist_ok=True)
        with baseline.open('xb') as stream:np.savez_compressed(stream,**values)
        result=dict(captured=True,arrays=len(values),all_finite=all(np.isfinite(v).all() for v in values.values()))
    else:
        with np.load(baseline) as saved:
            equal=set(saved.files)==set(values) and all(saved[k].shape==values[k].shape and
                np.array_equal(saved[k].view(np.uint64),values[k].view(np.uint64)) for k in saved.files)
        result=dict(captured=False,arrays=len(values),bitwise_equal=bool(equal))
    result.update(scope='default host API values; no solve or physical acceptance',native_images=images,
                  baseline=str(baseline.resolve()),baseline_sha256=hashlib.sha256(baseline.read_bytes()).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with receipt.open('x') as stream:stream.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
    return 0 if result.get('all_finite',result.get('bitwise_equal',False)) else 1


if __name__=='__main__':raise SystemExit(main())
