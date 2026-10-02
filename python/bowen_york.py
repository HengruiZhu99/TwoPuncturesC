"""Serial ctypes adapter for the BY equation system; explicit library paths only.

BY parameters are process-global. This adapter permits one live context per
process (even across backend instances); use fresh processes for parallel work.
P and S are lab-frame BY inputs, rather than Kerr seed velocities/rest spins.
"""
import ctypes as C
import hashlib
from pathlib import Path
import threading
import numpy as np
from native_loader import verify_image,validate_krylov
REAL_KEYS=frozenset(['par_b', 'par_m_plus', 'par_m_minus', 'target_M_plus', 'target_M_minus', 'par_P_plus1', 'par_P_plus2', 'par_P_plus3', 'par_P_minus1', 'par_P_minus2', 'par_P_minus3', 'par_S_plus1', 'par_S_plus2', 'par_S_plus3', 'par_S_minus1', 'par_S_minus2', 'par_S_minus3', 'center_offset1', 'center_offset2', 'center_offset3', 'Newton_tol', 'TP_linear_rtol', 'TP_epsilon', 'TP_Tiny', 'TP_Extend_Radius', 'adm_tol', 'initial_lapse_psi_exponent'])
INTEGER_KEYS=frozenset(['give_bare_mass', 'npoints_A', 'npoints_B', 'npoints_phi', 'TP_preconditioner', 'TP_linear_relative', 'TP_krylov_solver', 'TP_krylov_maxit', 'TP_krylov_restart', 'Newton_maxit', 'solve_momentum_constraint', 'use_external_initial_guess', 'do_residuum_debug_output', 'do_initial_debug_output', 'do_solution_file_output', 'do_bam_file_output', 'grid_setup_method', 'initial_lapse', 'conformal_state', 'swap_xz', 'multiply_old_lapse', 'verbose'])
PTR=C.POINTER(C.c_double)
def ptr(a):return a.ctypes.data_as(PTR)
class Derivs(C.Structure):
    _fields_=[('size',C.c_int)]+[(name,PTR) for name in ('d0','d1','d2','d3','d11','d12','d13','d22','d23','d33')]
class SolverStats(C.Structure):
    _fields_=[(key,C.c_int) for key in ('newton_iterations','krylov_iterations','jvp_applications','preconditioner_applications','relaxation_sweeps','modal_factorizations','linear_failures','modal_failures')]+[(key,C.c_double) for key in ('last_linear_target','last_true_linear_residual','last_relative_linear_residual')]
class InitialData(C.Structure):
    _fields_=[('F',PTR)]+[(name,C.POINTER(Derivs)) for name in ('u','v','cf_v')]+[('ntotal',C.c_int)]

class Backend:
    def __init__(self,library):
        path=Path(library)
        if not path.is_absolute():raise ValueError('native library path must be absolute')
        self.path=path.resolve(strict=True);self.lib=C.CDLL(str(self.path))
        self.loaded_sha256=verify_image(self.lib,'TwoPunctures_make_initial_data',self.path)
        signatures={'TwoPunctures_params_set_default':(None,[]),
          'TwoPunctures_params_reset':(C.c_int,[]),
          'TwoPunctures_params_set_Real':(None,[C.c_char_p,C.c_double]),
          'TwoPunctures_params_set_Int':(None,[C.c_char_p,C.c_int]),
          'TwoPunctures_make_initial_data':(C.c_void_p,[]),
          'TwoPunctures_finalise':(None,[C.c_void_p]),
          'TwoPunctures_diagnostics':(C.c_int,[C.c_void_p,PTR,PTR,PTR]),
          'TwoPunctures_sample_points':(C.c_int,[C.c_void_p,C.c_int,PTR,PTR,PTR,PTR,PTR]),
          'TP_solver_get_statistics':(None,[C.POINTER(SolverStats)])}
        for name,(ret,args) in signatures.items():
            if hasattr(self.lib,name):
                f=getattr(self.lib,name);f.restype=ret;f.argtypes=args
    def config(self):
        return dict(real=dict(par_b=3.,par_m_plus=.6,par_m_minus=.4,Newton_tol=1e-12),
                    integer=dict(npoints_A=40,npoints_B=80,npoints_phi=16,Newton_maxit=24,
                                 give_bare_mass=1,use_external_initial_guess=0,grid_setup_method=1,verbose=0))
    def create(self,config):return Solution(self,config)
    def library_sha256(self):
        if hashlib.sha256(self.path.read_bytes()).hexdigest()!=self.loaded_sha256:
            raise ValueError('native library file changed after loading; start a fresh process')
        return self.loaded_sha256

class Solution:
    _lock=threading.Lock()
    _poisoned=False
    def __init__(self,backend,config):
        self.context=None;self._owned=False;self._parameters_seeded=False
        if not self._lock.acquire(blocking=False):raise ValueError('BY permits one live context per process')
        self._owned=True;self.backend=backend
        try:
            if self._poisoned:raise ValueError('BY setup failed in an archived library; start a fresh process')
            if set(config)-{'real','integer'}:raise ValueError('BY config requires real/integer parameter dictionaries')
            self.config={kind:dict(config.get(kind,{})) for kind in ('real','integer')}
            if set(self.config['real'])-REAL_KEYS or set(self.config['integer'])-INTEGER_KEYS:
                raise ValueError('unknown BY parameter key or wrong parameter type')
        except Exception:self.close();raise
    def __enter__(self):return self
    def __exit__(self,*args):self.close()
    def __del__(self):self.close()
    def close(self):
        if self.context:self.backend.lib.TwoPunctures_finalise(self.context);self.context=None;self._parameters_seeded=False
        if self._parameters_seeded:
            if hasattr(self.backend.lib,'TwoPunctures_params_reset'):self.backend.lib.TwoPunctures_params_reset()
            else:Solution._poisoned=True
            self._parameters_seeded=False
        if self._owned:self._lock.release();self._owned=False
    def solve(self,linear_rtol=None,krylov=None,preconditioner='lines',max_krylov=100,restart=64):
        if not self._owned:raise ValueError('closed BY context')
        if self.context:raise ValueError('BY context already solved; create a fresh context')
        validate_krylov(krylov,linear_rtol)
        if preconditioner not in ('lines','modal'):raise ValueError('preconditioner must be lines or modal')
        if not isinstance(max_krylov,int) or not 1<=max_krylov<=100000:raise ValueError('invalid max_krylov')
        if not isinstance(restart,int) or not 1<=restart<=4096:raise ValueError('invalid restart')
        lib=self.backend.lib;new=hasattr(lib,'PK_solve')
        if not new and (krylov is not None or max_krylov!=100 or restart!=64):
            raise ValueError('library lacks selectable Krylov API')
        if not hasattr(lib,'TP_solver_get_statistics') and (linear_rtol is not None or preconditioner!='lines'):
            raise ValueError('library lacks solver options')
        reals=dict(self.config['real']);integers=dict(self.config['integer'])
        if set(reals)-REAL_KEYS or set(integers)-INTEGER_KEYS:raise ValueError('unknown BY parameter key or wrong parameter type')
        if hasattr(lib,'TP_solver_get_statistics'):
            integers.update(TP_preconditioner=int(preconditioner=='modal'),TP_linear_relative=int(linear_rtol is not None))
            reals['TP_linear_rtol']=.001 if linear_rtol is None else linear_rtol
        if new:integers.update(TP_krylov_solver=0 if krylov=='gmres' else 1,TP_krylov_maxit=max_krylov,TP_krylov_restart=restart)
        # Convert everything before mutating the global native table.
        reals={key.encode():float(value) for key,value in reals.items()}
        integers={key.encode():int(value) for key,value in integers.items()}
        if not all(np.isfinite(value) for value in reals.values()):raise ValueError('BY real parameters must be finite')
        if self._parameters_seeded and not hasattr(lib,'TwoPunctures_params_reset'):
            raise ValueError('failed archived BY setup; start a fresh process')
        if hasattr(lib,'TwoPunctures_params_reset') and lib.TwoPunctures_params_reset():
            raise ValueError('native BY data already live outside this adapter')
        self._parameters_seeded=True
        lib.TwoPunctures_params_set_default()
        for key,value in reals.items():lib.TwoPunctures_params_set_Real(key,value)
        for key,value in integers.items():lib.TwoPunctures_params_set_Int(key,value)
        self._resolved_config=dict(real={key.decode():value for key,value in reals.items()},integer={key.decode():value for key,value in integers.items()})
        self.context=lib.TwoPunctures_make_initial_data()
        if not self.context:raise ValueError('BY allocation/solve failed')
        self.resolved_options=dict(system='bowen_york',krylov=krylov or 'bicgstab',preconditioner=preconditioner,linear_rtol=linear_rtol,max_krylov=max_krylov,restart=restart)
        return self.diagnostics()
    def _check(self):
        if not self.context:raise ValueError('BY context not solved or closed')
    def diagnostics(self):
        self._check();residual=C.c_double();energy=C.c_double();masses=(C.c_double*2)()
        status=self.backend.lib.TwoPunctures_diagnostics(self.context,C.byref(residual),C.byref(energy),masses)
        tolerance=float(self._resolved_config['real'].get('Newton_tol',1e-10))
        target_ok=True
        if not self._resolved_config['integer'].get('give_bare_mass',1):
            targets=[float(self._resolved_config['real'].get('target_M_'+side,.5)) for side in ('plus','minus')]
            target_ok=max(abs(masses[i]-targets[i]) for i in range(2))<=float(self._resolved_config['real'].get('adm_tol',1e-10))
        converged=bool(status==0 and np.isfinite(residual.value) and residual.value<=tolerance and target_ok)
        return dict(status=0 if converged else 1,native_status=status,converged=converged,residual_linf=residual.value,adm_energy=energy.value,puncture_masses=list(masses),options=self.resolved_options.copy())
    def work_statistics(self):
        self._check()
        if not hasattr(self.backend.lib,'TP_solver_get_statistics'):return None
        s=SolverStats();self.backend.lib.TP_solver_get_statistics(C.byref(s))
        return {key:getattr(s,key) for key,_ in s._fields_}
    def unknowns(self):
        self._check();data=C.cast(self.context,C.POINTER(InitialData)).contents
        return np.ctypeslib.as_array(data.v.contents.d0,(data.ntotal,)).copy()
    def sample(self,xyz):
        self._check();points=np.ascontiguousarray(xyz,dtype=float).reshape(-1,3);n=len(points)
        lapse=np.empty(n);psi=np.empty(n);g6=np.empty((n,6));k6=np.empty_like(g6)
        if self.backend.lib.TwoPunctures_sample_points(self.context,n,ptr(points),ptr(lapse),ptr(psi),ptr(g6),ptr(k6)):
            raise ValueError('BY physical sampling failed')
        gamma=np.empty((n,3,3));K=np.empty_like(gamma)
        for q,(i,j) in enumerate(((0,0),(0,1),(0,2),(1,1),(1,2),(2,2))):
            gamma[:,i,j]=gamma[:,j,i]=g6[:,q];K[:,i,j]=K[:,j,i]=k6[:,q]
        return dict(gamma=gamma,Kij=K,psi=psi,lapse=lapse,attenuation=np.ones(n))
