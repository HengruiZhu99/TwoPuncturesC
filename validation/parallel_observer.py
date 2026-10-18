"""Ordered, read-only native sampling for the independent physical observer.

ctypes releases the GIL during native sampling. Warm the coefficient cache
before starting workers; all subsequent native work uses immutable arrays.
The same read-only sampler is already used by AthenaK's parallel geometry.
"""
from concurrent.futures import ThreadPoolExecutor
import numpy as np

class ParallelObserver:
    def __init__(self,solution,threads=1):
        if not isinstance(threads,int) or not 1<=threads<=32:raise ValueError('observer threads must lie in[1,32]')
        self.solution=solution;self.threads=threads;self.pool=None
    def __enter__(self):
        # Build coefficients serially. Zero points still calls the native
        # sampler but does not build its cache, so use one safe exterior point.
        centers=np.array([h.center for h in self.solution.config.hole]);origin=np.mean(centers,axis=0)
        radius=1+np.max(np.linalg.norm(centers-origin,axis=1))
        self.solution.sample([origin+radius*np.array([4.,3.,2.])])
        if self.threads>1:self.pool=ThreadPoolExecutor(max_workers=self.threads)
        return self
    def __exit__(self,*args):
        if self.pool:self.pool.shutdown(wait=True)
    def __call__(self,xyz):
        xyz=np.ascontiguousarray(xyz,dtype=float).reshape(-1,3)
        if self.pool is None or len(xyz)<self.threads:return self.solution.sample(xyz)
        pieces=list(self.pool.map(self.solution.sample,np.array_split(xyz,self.threads)))
        return {key:np.concatenate([part[key] for part in pieces]) for key in pieces[0]}
