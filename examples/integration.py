"""Isolated array export for later Lazarus integration; no project imports."""
import argparse,json
from pathlib import Path
import numpy as np
from hispid import Backend
from configs import moderate,as_dict

def main():
    p=argparse.ArgumentParser();p.add_argument('--library',required=True)
    p.add_argument('--n',type=int,default=40);p.add_argument('--nphi',type=int,default=20)
    p.add_argument('--output',default='validation/raw/integration.npz');a=p.parse_args()
    b=Backend(a.library);c=moderate(b,a.n,a.nphi)
    x=np.array([[3.5,.3,.2],[-2.5,-.3,.2],[8.,1.,-.7]])
    with b.create(c) as s:
        d=s.solve();fields=s.sample(x)
        output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(output,xyz=x,**fields)
        output.with_suffix('.json').write_text(json.dumps(dict(config=as_dict(c),diagnostics=d,library=str(b.path),physical_acceptance='run independent validation before use'),indent=2)+'\n')
        print(json.dumps(d,indent=2))
    return 0 if d['status']==0 else 1

if __name__=='__main__':raise SystemExit(main())
