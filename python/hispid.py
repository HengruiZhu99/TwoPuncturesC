"""Explicit-path ctypes adapter to the isolated HiSpID native backend.

No library discovery, global solver parameters, or Lazarus installation edits.
Requires NumPy. Use ``Backend('/absolute/path/build-hispid/libHiSpID.so')``.
"""
from __future__ import annotations
import ctypes as C
import hashlib
import sys
from pathlib import Path
import numpy as np

D3=C.c_double*3
D9=C.c_double*9

class Hole(C.Structure):
    _fields_=[('mass',C.c_double),('center',D3),('spin',D3),('velocity',D3)]
    def __init__(self,mass=0,center=(0,0,0),spin=(0,0,0),velocity=(0,0,0)):
        super().__init__(mass,D3(*center),D3(*spin),D3(*velocity))

class Config(C.Structure):
    _fields_=[('hole',Hole*2),('n',C.c_int*3),('conformal_choice',C.c_int),
              ('inner_flatten',C.c_int),('omega',C.c_double*2),
              ('attenuation_power',C.c_int),('inner_min',C.c_double*2),
              ('inner_max',C.c_double*2),('far_radius',C.c_double),
              ('tolerance',C.c_double),('max_newton',C.c_int),
              ('max_krylov',C.c_int),('krylov_restart',C.c_int),('memory_limit_mib',C.c_int)]

class Point(C.Structure):
    _fields_=[('gamma',D9),('Kij',D9),('psi',C.c_double),('conformal_metric',D9),
              ('Atilde',D9),('mean_curvature',C.c_double),
              ('correction',C.c_double*4),('attenuation',C.c_double)]

class Diagnostics(C.Structure):
    _fields_=[('converged',C.c_int),('newton_iterations',C.c_int),
              ('krylov_iterations',C.c_int),('npoints',C.c_int),
              ('scaled_linf',C.c_double*4),('unscaled_linf',C.c_double*4),
              ('seconds',C.c_double)]
    def as_dict(self):
        return {k:list(getattr(self,k)) if isinstance(getattr(self,k),C.Array)
                else getattr(self,k) for k,_ in self._fields_}

PTR=C.POINTER(C.c_double)
def ptr(a):return a.ctypes.data_as(PTR)
def unpack(out):
    return {name:np.array([list(getattr(p,name)) if isinstance(getattr(p,name),C.Array)
                          else getattr(p,name) for p in out])
            for name,_ in Point._fields_}

class Backend:
    _loaded_builds = {}
    def __init__(self,library):
        path=Path(library)
        if not path.is_absolute():raise ValueError('native library path must be absolute')
        self.path=path.resolve(strict=True)
        self.lib=C.CDLL(str(self.path))
        # dyld identifies copied libraries by their embedded install name.
        # Loading two archived builds in one process can silently reuse the
        # first image. Check the function's actual image before any API call.
        address=C.cast(self.lib.HiSpID_default_config,C.c_void_p).value
        if sys.platform=='darwin' or sys.platform.startswith('linux'):
            class DlInfo(C.Structure):
                _fields_=[('filename',C.c_char_p),('base',C.c_void_p),
                          ('symbol',C.c_char_p),('symbol_address',C.c_void_p)]
            dladdr=C.CDLL(None).dladdr
            dladdr.argtypes=[C.c_void_p,C.POINTER(DlInfo)];dladdr.restype=C.c_int
            info=DlInfo()
            if not dladdr(address,C.byref(info)) or not info.filename:
                raise ValueError('cannot verify loaded native library image')
            actual=Path(info.filename.decode()).resolve(strict=True)
            if actual!=self.path:
                raise ValueError(f'native loader reused {actual}; compare distinct builds in separate processes')
        digest=hashlib.sha256(self.path.read_bytes()).hexdigest()
        previous=self._loaded_builds.get(address)
        if previous is not None and previous!=(self.path,digest):
            raise ValueError('loaded native library differs from on-disk build; start a fresh process')
        self._loaded_builds[address]=(self.path,digest)
        self.loaded_sha256=digest
        api={'HiSpID_default_config':(None,[C.POINTER(Config)]),
             'HiSpID_create':(C.c_void_p,[C.POINTER(Config)]),
             'HiSpID_solve':(C.c_int,[C.c_void_p]),
             'HiSpID_destroy':(None,[C.c_void_p]),
             'HiSpID_seed':(C.c_int,[C.POINTER(Hole),C.c_int,C.c_int,PTR,C.POINTER(Point)]),
             'HiSpID_sample':(C.c_int,[C.c_void_p,C.c_int,PTR,C.POINTER(Point)]),
             'HiSpID_diagnostics':(C.c_int,[C.c_void_p,C.POINTER(Diagnostics)]),
             'HiSpID_get_unknowns':(C.c_int,[C.c_void_p,PTR,C.c_int]),
             'HiSpID_set_unknowns':(C.c_int,[C.c_void_p,PTR,C.c_int]),
             'HiSpID_residual':(C.c_int,[C.c_void_p,PTR,PTR]),
             'HiSpID_jvp':(C.c_int,[C.c_void_p,PTR,PTR,PTR]),
             'HiSpID_equation_samples':(C.c_int,[C.c_void_p,PTR,PTR,PTR,PTR]),
             'HiSpID_operators':(C.c_int,[C.POINTER(Config),PTR,PTR,PTR]),
             'HiSpID_charges':(C.c_int,[C.c_void_p,PTR,C.c_double,C.c_int,C.c_int,PTR]),
             'HiSpID_last_error':(C.c_char_p,[])}
        for name,(ret,args) in api.items():
            f=getattr(self.lib,name);f.restype=ret;f.argtypes=args
        # Archived libraries remain loadable for explicit API migration checks.
        optional={'HiSpID_create_sampler':(C.c_void_p,[C.POINTER(Config)]),
                  'HiSpID_sample_with_derivatives':(C.c_int,[C.c_void_p,C.c_int,PTR,C.POINTER(Point),PTR])}
        for name,(ret,args) in optional.items():
            if hasattr(self.lib,name):
                f=getattr(self.lib,name);f.restype=ret;f.argtypes=args
    def error(self):return self.lib.HiSpID_last_error().decode()
    def library_sha256(self):
        if hashlib.sha256(self.path.read_bytes()).hexdigest()!=self.loaded_sha256:
            raise ValueError('native library file changed after loading; start a fresh process')
        return self.loaded_sha256
    def config(self):
        c=Config();self.lib.HiSpID_default_config(C.byref(c));return c
    def seed(self,hole,xyz,choice=1):
        x=np.ascontiguousarray(xyz,dtype=float).reshape(-1,3);out=(Point*len(x))()
        r=self.lib.HiSpID_seed(C.byref(hole),choice,len(x),ptr(x),out)
        if r:raise ValueError(self.error() or 'invalid seed input')
        return unpack(out)
    def operators(self,config,xyz,jets):
        x=np.ascontiguousarray(xyz,dtype=float).reshape(3);j=np.ascontiguousarray(jets,dtype=float).reshape(40);out=np.empty(5)
        if self.lib.HiSpID_operators(C.byref(config),ptr(x),ptr(j),ptr(out)):raise ValueError(self.error())
        return out
    def create(self,config):return Solution(self,config)
    def create_sampler(self,config):return Solution(self,config,sampler_only=True)

class Solution:
    def __init__(self,backend,config,sampler_only=False):
        self.backend=backend;self.config=Config.from_buffer_copy(config)
        self.context=None
        name='HiSpID_create_sampler' if sampler_only else 'HiSpID_create'
        if not hasattr(backend.lib,name):raise ValueError('library does not support sampling-only contexts')
        self.context=getattr(backend.lib,name)(C.byref(config))
        if not self.context:raise ValueError(backend.error())
        self.size=4*int(np.prod(list(config.n)))
    def __enter__(self):return self
    def __exit__(self,*args):self.close()
    def close(self):
        if self.context:self.backend.lib.HiSpID_destroy(self.context);self.context=None
    def __del__(self):self.close()
    def _check(self):
        if not self.context:raise ValueError('closed HiSpID context')
    def solve(self):
        self._check();r=self.backend.lib.HiSpID_solve(self.context)
        d=self.diagnostics();d['status']=r;d['error']=self.backend.error() if r else '';return d
    def diagnostics(self):
        self._check();d=Diagnostics();self.backend.lib.HiSpID_diagnostics(self.context,C.byref(d));return d.as_dict()
    def sample(self,xyz):
        self._check();x=np.ascontiguousarray(xyz,dtype=float).reshape(-1,3);out=(Point*len(x))()
        if self.backend.lib.HiSpID_sample(self.context,len(x),ptr(x),out):raise ValueError(self.backend.error())
        return unpack(out)
    def sample_with_derivatives(self,xyz):
        self._check()
        if not hasattr(self.backend.lib,'HiSpID_sample_with_derivatives'):
            raise ValueError('library does not support physical metric gradients')
        x=np.ascontiguousarray(xyz,dtype=float).reshape(-1,3);out=(Point*len(x))();dg=np.empty((len(x),3,3,3))
        if self.backend.lib.HiSpID_sample_with_derivatives(self.context,len(x),ptr(x),out,ptr(dg)):
            raise ValueError(self.backend.error())
        values=unpack(out);values['dgamma']=dg;return values
    def unknowns(self):
        self._check();out=np.empty(self.size);self.backend.lib.HiSpID_get_unknowns(self.context,ptr(out),self.size);return out
    def set_unknowns(self,values):
        self._check();v=np.ascontiguousarray(values,dtype=float)
        if self.backend.lib.HiSpID_set_unknowns(self.context,ptr(v),v.size):raise ValueError('invalid unknown vector')
    def residual(self,values):
        self._check();v=np.ascontiguousarray(values,dtype=float)
        if v.size!=self.size:raise ValueError('wrong unknown count')
        out=np.empty(self.size)
        if self.backend.lib.HiSpID_residual(self.context,ptr(v),ptr(out)):raise ValueError(self.backend.error() or 'residual evaluation failed')
        return out
    def jvp(self,base,direction):
        self._check();b=np.ascontiguousarray(base,dtype=float);d=np.ascontiguousarray(direction,dtype=float)
        if b.size!=self.size or d.size!=self.size:raise ValueError('wrong unknown count')
        out=np.empty(self.size)
        if self.backend.lib.HiSpID_jvp(self.context,ptr(b),ptr(d),ptr(out)):raise ValueError(self.backend.error() or 'Jacobian evaluation failed')
        return out
    def charges(self,radius,center=(0,0,0),ntheta=16,nphi=32):
        self._check();x=np.ascontiguousarray(center,dtype=float).reshape(3);out=np.empty(7)
        if self.backend.lib.HiSpID_charges(self.context,ptr(x),radius,ntheta,nphi,ptr(out)):raise ValueError(self.backend.error())
        return out
    def equation_samples(self):
        self._check();n=self.size//4;x=np.empty((n,3));g=np.empty(n);psi=np.empty(n);HM=np.empty((n,4))
        if self.backend.lib.HiSpID_equation_samples(self.context,ptr(x),ptr(g),ptr(psi),ptr(HM)):raise ValueError(self.backend.error())
        return dict(xyz=x,attenuation=g,psi=psi,physical_equivalent=HM)
