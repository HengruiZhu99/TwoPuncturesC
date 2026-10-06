"""Explicit-path ctypes adapter to the isolated HiSpID native backend.

No library discovery, global solver parameters, or Lazarus installation edits.
Requires NumPy. Use ``Backend('/absolute/path/build-hispid/libHiSpID.so')``.
"""
from __future__ import annotations
import ctypes as C
import hashlib
from pathlib import Path
import numpy as np
from native_loader import verify_image,verify_puncture_dependencies,validate_krylov

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

class SolveOptions(C.Structure):
    _fields_=[('struct_size',C.c_int),('krylov',C.c_int),('linear_rtol',C.c_double)]

class SetupStatistics(C.Structure):
    _fields_=[('struct_size',C.c_int),('geometry_execution',C.c_int),('scalar_digits',C.c_int),
              ('spectral_seconds',C.c_double),('geometry_seconds',C.c_double),('coefficient_seconds',C.c_double)]

PTR=C.POINTER(C.c_double)
def ptr(a):return a.ctypes.data_as(PTR)
def unpack(out):
    return {name:np.array([list(getattr(p,name)) if isinstance(getattr(p,name),C.Array)
                          else getattr(p,name) for p in out])
            for name,_ in Point._fields_}

SEED_FAMILIES={'qi':0,'trumpet_r0_m':1}
def seed_family_code(name):
    if name not in SEED_FAMILIES:raise ValueError('unsupported seed family: '+str(name))
    return SEED_FAMILIES[name]

class Backend:
    def __init__(self,library):
        path=Path(library)
        if not path.is_absolute():raise ValueError('native library path must be absolute')
        self.path=path.resolve(strict=True)
        self.lib=C.CDLL(str(self.path))
        self.loaded_sha256=verify_image(self.lib,'HiSpID_default_config',self.path)
        self.dependency_images=verify_puncture_dependencies(self.lib,self.path)
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
        optional={'HiSpID_create_with_seed_family':(C.c_void_p,[C.POINTER(Config),C.c_int,C.c_int,C.c_int,C.c_int]),
                  'HiSpID_seed_family':(C.c_int,[C.c_void_p]),
                  'HiSpID_seed_with_family':(C.c_int,[C.POINTER(Hole),C.c_int,C.c_int,PTR,C.POINTER(Point),C.c_int,C.c_int]),
                  'HiSpID_operators_with_seed_family':(C.c_int,[C.POINTER(Config),PTR,PTR,PTR,C.c_int]),
                  'HiSpID_work_statistics':(C.c_int,[C.c_void_p,C.POINTER(C.c_int)]),
                  'HiSpID_create_with_execution':(C.c_void_p,[C.POINTER(Config),C.c_int]),
                  'HiSpID_create_with_geometry':(C.c_void_p,[C.POINTER(Config),C.c_int,C.c_int]),
                  'HiSpID_seed_with_execution':(C.c_int,[C.POINTER(Hole),C.c_int,C.c_int,PTR,C.POINTER(Point),C.c_int]),
                  'HiSpID_setup_statistics':(C.c_int,[C.c_void_p,C.POINTER(SetupStatistics)]),
                  'HiSpID_default_solve_options':(None,[C.POINTER(SolveOptions)]),
                  'HiSpID_solve_with_options':(C.c_int,[C.c_void_p,C.POINTER(SolveOptions)]),
                  'HiSpID_set_axisymmetric':(C.c_int,[C.c_void_p,C.c_int]),
                  'HiSpID_resolved_solve_options':(C.c_int,[C.c_void_p,C.POINTER(SolveOptions)]),
                  'HiSpID_linear_history':(C.c_int,[C.c_void_p,C.c_int,PTR]),
                  'HiSpID_solve_with_forcing':(C.c_int,[C.c_void_p,C.c_double]),
                  'HiSpID_unknown_parameterization':(C.c_char_p,[]),
                  'HiSpID_residual_scaling':(C.c_char_p,[]),
                  'HiSpID_collocation_maps':(C.c_int,[PTR]),
                  'HiSpID_create_sampler':(C.c_void_p,[C.POINTER(Config)]),
                  'HiSpID_sample_with_derivatives':(C.c_int,[C.c_void_p,C.c_int,PTR,C.POINTER(Point),PTR])}
        for name,(ret,args) in optional.items():
            if hasattr(self.lib,name):
                f=getattr(self.lib,name);f.restype=ret;f.argtypes=args
    def parameterization(self):
        if hasattr(self.lib,'HiSpID_unknown_parameterization'):
            return self.lib.HiSpID_unknown_parameterization().decode('ascii')
        return 'W_plus_Aminus1_V'
    def residual_scaling(self):
        return self.lib.HiSpID_residual_scaling().decode('ascii') if hasattr(self.lib,'HiSpID_residual_scaling') else 'sin6_alpha_beta'
    def parameterization_description(self):
        if self.parameterization().startswith('modal_P_C2tauC4_map_v5_'):
            return 'C2 modal P with homogeneous axis endpoint constraints for m>=5 (C4 regularity)'
        if self.parameterization().startswith('modal_P_C4prolate_map_v4_'):
            return 'modal P: u=W-2(1-t)q^r P, t=a^2, q=a sin(R), parity cap r<=6'
        value='u=W+(A-1)V, W=sum((1-F)*(psi_seed-1))'
        if self.parameterization() in ('modal_P_C2prolate_v1','modal_P_C2prolate_mapped_v2') or self.parameterization().startswith('modal_P_C2prolate_map_v3_'):
            return 'modal P: u=W-2(1-t)q^r P, t=a^2, q=a sin(R), parity cap r<=4'
        return value+'; Cartesian C2 axis cardinal basis' if self.parameterization().endswith('_C2axis') else value
    def parameterization_maps(self):
        if hasattr(self.lib,'HiSpID_collocation_maps'):
            values=np.empty(2)
            if self.lib.HiSpID_collocation_maps(ptr(values)):raise ValueError('native collocation-map query failed')
            return dict(radial_stretch=float(values[0]),angular_stretch=float(values[1]))
        if self.parameterization()=='modal_P_C2prolate_mapped_v2':
            return dict(radial_stretch=.2,angular_stretch=2.)
        if self.parameterization()=='modal_P_C2prolate_v1':
            return dict(radial_stretch=1.,angular_stretch=0.)
        return None
    def error(self):return self.lib.HiSpID_last_error().decode()
    def library_sha256(self):
        if hashlib.sha256(self.path.read_bytes()).hexdigest()!=self.loaded_sha256:
            raise ValueError('native library file changed after loading; start a fresh process')
        for path,digest in self.dependency_images.items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:raise ValueError('native dependency changed after loading; start a fresh process')
        return self.loaded_sha256
    def config(self):
        c=Config();self.lib.HiSpID_default_config(C.byref(c));return c
    def seed(self,hole,xyz,choice=1,execution='reference',seed_family='qi'):
        x=np.ascontiguousarray(xyz,dtype=float).reshape(-1,3);out=(Point*len(x))()
        from execution import select
        code=select(self.lib,execution)
        family=seed_family_code(seed_family)
        if family:
            if not hasattr(self.lib,'HiSpID_seed_with_family'):raise ValueError('library lacks seed-family support')
            r=self.lib.HiSpID_seed_with_family(C.byref(hole),choice,len(x),ptr(x),out,code,family)
        elif code:
            if not hasattr(self.lib,'HiSpID_seed_with_execution'):raise ValueError('library lacks execution-space seed export')
            r=self.lib.HiSpID_seed_with_execution(C.byref(hole),choice,len(x),ptr(x),out,code)
        else:r=self.lib.HiSpID_seed(C.byref(hole),choice,len(x),ptr(x),out)
        if r:raise ValueError(self.error() or 'invalid seed input')
        return unpack(out)
    def operators(self,config,xyz,jets,seed_family=None):
        x=np.ascontiguousarray(xyz,dtype=float).reshape(3);j=np.ascontiguousarray(jets,dtype=float).reshape(40);out=np.empty(5)
        family=seed_family_code(getattr(config,'seed_family','qi') if seed_family is None else seed_family)
        if family:
            if not hasattr(self.lib,'HiSpID_operators_with_seed_family'):raise ValueError('library lacks seed-family operators')
            code=self.lib.HiSpID_operators_with_seed_family(C.byref(config),ptr(x),ptr(j),ptr(out),family)
        else:code=self.lib.HiSpID_operators(C.byref(config),ptr(x),ptr(j),ptr(out))
        if code:raise ValueError(self.error())
        return out
    def create(self,config,execution='reference',geometry='host',seed_family=None):return Solution(self,config,execution=execution,geometry=geometry,seed_family=seed_family)
    def create_sampler(self,config,seed_family=None):return Solution(self,config,sampler_only=True,seed_family=seed_family)

class Solution:
    def __init__(self,backend,config,sampler_only=False,execution='reference',geometry='host',seed_family=None):
        self.backend=backend;self.config=Config.from_buffer_copy(config)
        self.context=None
        configured=getattr(config,'seed_family','qi')
        if seed_family is not None and hasattr(config,'seed_family') and seed_family!=configured:
            raise ValueError('seed-family override conflicts with configuration metadata')
        self.seed_family=configured if seed_family is None else seed_family
        family=seed_family_code(self.seed_family);self.config.seed_family=self.seed_family
        name='HiSpID_create_sampler' if sampler_only else 'HiSpID_create'
        if not hasattr(backend.lib,name):raise ValueError('library does not support sampling-only contexts')
        from execution import select
        code=select(backend.lib,execution)
        if geometry not in ('host','execution'):raise ValueError('geometry must be host or execution')
        if sampler_only and (code or geometry!='host'):raise ValueError('sampling-only contexts use host geometry and CPU execution')
        if family:
            if not hasattr(backend.lib,'HiSpID_create_with_seed_family'):raise ValueError('library lacks seed-family support')
            self.context=backend.lib.HiSpID_create_with_seed_family(C.byref(config),family,code,int(geometry=='execution'),int(sampler_only))
        elif geometry=='execution':
            if not code or sampler_only:raise ValueError('execution geometry requires a Kokkos solving context')
            if not hasattr(backend.lib,'HiSpID_create_with_geometry'):raise ValueError('library lacks execution-space geometry setup')
            self.context=backend.lib.HiSpID_create_with_geometry(C.byref(config),code,1)
        elif code:
            if sampler_only:raise ValueError('sampling-only contexts use CPU execution')
            if not hasattr(backend.lib,'HiSpID_create_with_execution'):raise ValueError('library lacks execution-aware contexts')
            self.context=backend.lib.HiSpID_create_with_execution(C.byref(config),code)
        else:self.context=getattr(backend.lib,name)(C.byref(config))
        self.execution=execution;self.geometry=geometry
        if not self.context:raise ValueError(backend.error())
        self.size=4*int(np.prod(list(config.n)))
    def __enter__(self):return self
    def __exit__(self,*args):self.close()
    def close(self):
        if self.context:self.backend.lib.HiSpID_destroy(self.context);self.context=None
    def __del__(self):self.close()
    def _check(self):
        if not self.context:raise ValueError('closed HiSpID context')
    def setup_statistics(self):
        self._check()
        if not hasattr(self.backend.lib,'HiSpID_setup_statistics'):return None
        out=SetupStatistics();out.struct_size=C.sizeof(out)
        if self.backend.lib.HiSpID_setup_statistics(self.context,C.byref(out)):raise ValueError('native setup-statistics query failed')
        return {name:getattr(out,name) for name,_ in out._fields_ if name!='struct_size'}
    def solve(self,linear_rtol=None,krylov=None,axisymmetric=False):
        self._check()
        validate_krylov(krylov,linear_rtol,allow_lgmres=True)
        if not isinstance(axisymmetric,bool):raise ValueError('axisymmetric must be boolean')
        if hasattr(self.backend.lib,'HiSpID_set_axisymmetric'):
            if self.backend.lib.HiSpID_set_axisymmetric(self.context,int(axisymmetric)):
                raise ValueError(self.backend.error())
        elif axisymmetric:raise ValueError('library lacks axisymmetric solve API')
        if krylov is not None:
            if not hasattr(self.backend.lib,'HiSpID_solve_with_options'):raise ValueError('library lacks selectable Krylov API')
            options=SolveOptions(C.sizeof(SolveOptions),{'gmres':0,'bicgstab':1,'lgmres':2}[krylov],0 if linear_rtol is None else float(linear_rtol))
            r=self.backend.lib.HiSpID_solve_with_options(self.context,C.byref(options))
        elif linear_rtol is None:r=self.backend.lib.HiSpID_solve(self.context)
        else:
            if not hasattr(self.backend.lib,'HiSpID_solve_with_forcing'):raise ValueError('library lacks fixed forcing API')
            r=self.backend.lib.HiSpID_solve_with_forcing(self.context,float(linear_rtol))
        self.resolved_options=dict(system='hispid',krylov=krylov or 'gmres',linear_rtol=linear_rtol,preconditioner='modal',execution=self.execution,geometry=self.geometry,seed_family=self.seed_family)
        if axisymmetric:self.resolved_options['axisymmetric']=True
        self.resolved_options['native_verified']=False
        if hasattr(self.backend.lib,'HiSpID_resolved_solve_options'):
            actual=SolveOptions(C.sizeof(SolveOptions),0,0)
            if self.backend.lib.HiSpID_resolved_solve_options(self.context,C.byref(actual))==0:
                self.resolved_options.update(krylov={0:'gmres',1:'bicgstab',2:'lgmres'}[actual.krylov],linear_rtol=actual.linear_rtol or None,native_verified=True)
        d=self.diagnostics();d['status']=r;d['error']=self.backend.error() if r else '';return d
    def work_statistics(self):
        self._check()
        if not hasattr(self.backend.lib,'HiSpID_work_statistics'):return None
        values=(C.c_int*2)()
        if self.backend.lib.HiSpID_work_statistics(self.context,values):raise ValueError('work statistics unavailable')
        return dict(jvp_applications=values[0],preconditioner_applications=values[1])
    def linear_history(self):
        self._check()
        if not hasattr(self.backend.lib,'HiSpID_linear_history'):return None
        size=self.backend.lib.HiSpID_linear_history(self.context,0,None)
        if size<0:raise ValueError('linear history unavailable')
        out=np.empty((size,4))
        if size and self.backend.lib.HiSpID_linear_history(self.context,size,ptr(out))<0:raise ValueError('linear history unavailable')
        return out.tolist()
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
