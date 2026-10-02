"""Fresh-process original standalone TwoPunctures reference state worker.

Compatible with upstream ec563aeb, without using fork-only sampling/diagnostic
symbols. The original Cartesian scalar interpolator is evaluated identically.
"""
import argparse
import ctypes as C
import json
from pathlib import Path
import resource
import sys
import numpy as np
from benchmark_bowen_york import InitialData, save_by_state, verify_image, digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library',required=True);p.add_argument('--parameters',required=True)
    p.add_argument('--state',required=True);p.add_argument('--output',required=True)
    p.add_argument('--target-mass',action='store_true')
    args=p.parse_args();path=Path(args.library).resolve(strict=True);expected_sha=digest(path);lib=C.CDLL(str(path));verify_image(lib,'TwoPunctures_make_initial_data',path)
    record=json.loads(Path(args.parameters).read_text());reals=record['by_real_parameters'].copy();integers=record['by_integer_parameters'].copy()
    integers['verbose']=1
    if args.target_mass:
        integers['give_bare_mass']=0
        reals.update(target_M_plus=.6,target_M_minus=.4,adm_tol=1e-10)
    for name,restype,argtypes in [('TwoPunctures_params_set_default',None,[]),('TwoPunctures_params_set_Real',None,[C.c_char_p,C.c_double]),('TwoPunctures_params_set_Int',None,[C.c_char_p,C.c_int]),('TwoPunctures_make_initial_data',C.c_void_p,[])]:
        f=getattr(lib,name);f.restype=restype;f.argtypes=argtypes
    lib.TwoPunctures_params_set_default()
    for k,v in reals.items():lib.TwoPunctures_params_set_Real(k.encode(),v)
    for k,v in integers.items():lib.TwoPunctures_params_set_Int(k.encode(),v)
    pointer=lib.TwoPunctures_make_initial_data()
    if not pointer:raise RuntimeError('original solve failed')
    data=C.cast(pointer,C.POINTER(InitialData)).contents
    residual=float(np.max(np.abs(np.ctypeslib.as_array(data.F,(data.ntotal,)))))
    # Avoid fork-only sample symbols even when the worker is applied to a fork.
    save_by_state(lib,pointer,record['config']['n'],args.state,include_physical=False)
    # A fresh process owns/frees its address space on exit; upstream finalise has
    # different parameter ownership semantics, immaterial to this first solve.
    lib.params_get_real.argtypes=[C.c_char_p];lib.params_get_real.restype=C.c_double
    masses=[lib.params_get_real(name) for name in (b'par_m_plus',b'par_m_minus')]
    lib.PunctIntPolAtArbitPosition.argtypes=[C.c_int]*5+[C.POINTER(type(data.u.contents))]+[C.c_double]*3
    lib.PunctIntPolAtArbitPosition.restype=C.c_double
    grid=record['config']['n'];separation=reals['par_b']
    up=lib.PunctIntPolAtArbitPosition(0,1,*grid,data.v,separation,0,0)
    um=lib.PunctIntPolAtArbitPosition(0,1,*grid,data.v,-separation,0,0)
    end_masses=[(1+up)*masses[0]+masses[0]*masses[1]/(4*separation),
                (1+um)*masses[1]+masses[0]*masses[1]/(4*separation)]
    if digest(path)!=expected_sha:raise RuntimeError('reference library changed during solve')
    result=dict(final_bare_masses=masses,internal_end_adm_masses=end_masses,
                target_mass_errors=[abs(end_masses[0]-.6),abs(end_masses[1]-.4)] if args.target_mass else None,
                library_sha256=expected_sha,state_sha256=digest(args.state),residual_linf=residual,
                converged=bool(np.isfinite(residual) and residual<=reals['Newton_tol']),
                real_parameters=reals,integer_parameters=integers,loaded_image_verified=True)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print('REFERENCE',json.dumps(result),flush=True)
if __name__=='__main__':main()
