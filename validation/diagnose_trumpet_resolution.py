"""Directional residual replay of one retained polynomial; no new solves.

Three independent coordinate directions detect aliasing left by the coarse
collocation grid. Native equivalent constraints are diagnostics, not a physical
acceptance observer. Every replay retains the same continuous coefficients.
"""
import argparse,json,sys,hashlib,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Backend,Config
from checkpoint_export import read_checkpoint
from prolong import for_backend

def run(library,checkpoint,output):
    if output.exists():raise FileExistsError(output)
    b=Backend(str(library.resolve()));c,u,meta=read_checkpoint(checkpoint)
    if meta['source_library_sha256']!=b.loaded_sha256:raise ValueError('checkpoint image mismatch')
    result=dict(checkpoint=meta,driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),rows=[],binary_acceptance=False)
    old=list(c.n)
    for axis in range(3):
        shape=old.copy();shape[axis]*=2;fine=Config.from_buffer_copy(c);fine.seed_family=c.seed_family;fine.n[:]=shape
        values=for_backend(b,u,old,shape);start=time.monotonic()
        with b.create(fine,execution='kokkos',geometry='host') as s:
            s.set_unknowns(values);eq=s.equation_samples();F=eq['physical_equivalent'];x=eq['xyz'];g=eq['attenuation']
            distance=np.min([np.linalg.norm(x-np.array(h.center),axis=1)/h.mass for h in c.hole],axis=0)
            rows={}
            for label,mask in [('near',(distance<2)&(g==1)),('bulk',(distance>=2)&(np.linalg.norm(x,axis=1)<30)&(g==1))]:
                rows[label]=dict(count=int(mask.sum()),component_rms=np.sqrt(np.mean(F[mask]**2,axis=0)).tolist(),component_linf=np.max(abs(F[mask]),axis=0).tolist())
            result['rows'].append(dict(refined_axis=axis,resolution=shape,seconds=time.monotonic()-start,native_equivalent=rows))
        output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['rows'][-1]),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--library',type=Path,required=True);ap.add_argument('--checkpoint',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();run(a.library,a.checkpoint,a.output)
