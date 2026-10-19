import sys,json
from pathlib import Path
import numpy as np
root=Path.cwd();sys.path[:0]=[str(root/'python'),str(root/'validation')]
from hispid import Backend,Wave22
from kerr_wave_case import config
b=Backend(str(Path(sys.argv[1]).resolve()));c=config(b,8,.9999);c.hole[1].center[2]=-float(sys.argv[2])
w=Wave22(**json.loads((root/'tmp/new-map-polar-selection.json').read_text())['wave22'])
r=np.unique(np.r_[np.geomspace(.08,50,97),.1,.2,.25353544551647955,.2939671]);mu=np.array([-.95,-.6,0,.6,.95]);phi=np.pi/8
rays=np.stack([np.sqrt(1-mu*mu)*np.cos(phi),np.sqrt(1-mu*mu)*np.sin(phi),mu],axis=-1);x=(r[:,None,None]*rays[None,:,:]).reshape(-1,3)
with b.create_sampler(c,wave=w) as s:q=s.sample_with_derivatives(x)
np.savez_compressed(sys.argv[3],xyz=x,**{k:q[k] for k in ('gamma','Kij','dgamma')})
print(json.dumps(dict(library_sha256=b.library_sha256(),maps=b.parameterization_maps(),points=len(x))))
