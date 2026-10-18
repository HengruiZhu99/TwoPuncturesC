"""Serial/threaded field samples must be bit for bit identical and ordered."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation'),str(ROOT/'examples')]
from hispid import Backend
from trumpet_configs import qi_moderate_wide
from run_validation import points
from parallel_observer import ParallelObserver
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--threads',type=int,default=2);a=p.parse_args();b=Backend(str(a.library.resolve()));c=qi_moderate_wide(b,12,16)
with b.create_sampler(c) as s:
 n=4*np.prod(c.n);s.set_unknowns(1e-8*np.cos(.23*np.arange(n)));x,*_=points(c,True)
 old=s.sample(x)
 with ParallelObserver(s,a.threads) as observer:new=observer(x)
 assert old.keys()==new.keys()
 assert all(np.array_equal(old[k],new[k]) for k in old)
a.output.write_text(json.dumps(dict(passed=True,bitwise_equal=True,threads=a.threads,point_count=len(x),library_sha256=b.library_sha256()),indent=2)+'\n');print('all fields bitwise equal')
