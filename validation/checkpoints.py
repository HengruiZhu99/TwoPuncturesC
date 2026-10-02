"""Strict replay of checkpoints bound to their configuration and native SHA."""
import ctypes as C
import hashlib
import json
from pathlib import Path
import numpy as np
from hispid import Config,Hole

ROOT=Path(__file__).resolve().parents[1]
LEGACY_PARAMETERIZATION='u=W+(A-1)V, W=sum((1-F)*(psi_seed-1))'
LIFT_PARAMETERIZATION=LEGACY_PARAMETERIZATION+'; Cartesian C2 axis cardinal basis'
PARAMETERIZATION='modal P: u=W-2(1-t)q^r P, t=a^2, q=a sin(R), parity cap r<=4'

def library_sha(backend):
    return backend.library_sha256()

def select_record(case,resolution=None,nphi=None):
    result=json.loads((ROOT/'validation/results.json').read_text())[case]
    if resolution is None:return result['records'][-1]
    candidates=[r for r in result['records'] if r['resolution'][0]==resolution
                and (nphi is None or r['resolution'][2]==nphi)]
    if len(candidates)!=1:raise ValueError('select a unique checkpoint with --resolution and --nphi')
    return candidates[0]

def restore(backend,record):
    if record.get('library_sha256')!=library_sha(backend):
        raise ValueError('checkpoint/library mismatch; do not reinterpret old coefficients')
    token=backend.parameterization()
    if token!='W_plus_Aminus1_V' and record.get('unknown_parameterization_id')!=token:
        raise ValueError('checkpoint/native exact basis identifier mismatch')
    if record.get('collocation_maps')!=backend.parameterization_maps():
        raise ValueError('checkpoint/native collocation maps mismatch')
    if record.get('unknown_parameterization')!=backend.parameterization_description():
        raise ValueError('checkpoint/native continuous basis mismatch')
    return restore_payload(record,backend.config())

def restore_payload(record,cfg):
    """Decode typed payload only; caller must enforce its source SHA explicitly.

    This is used by isolated migration workers, never by normal replay.
    """
    if record.get('unknown_parameterization') not in (LEGACY_PARAMETERIZATION,LIFT_PARAMETERIZATION,PARAMETERIZATION):
        raise ValueError('unsupported checkpoint parameterization')
    for name,_ in Config._fields_:
        if name not in record['config']:
            if name=='memory_limit_mib':continue
            raise ValueError('missing checkpoint configuration field: '+name)
        value=record['config'][name]
        if name=='hole':
            for i,hole in enumerate(value):cfg.hole[i]=Hole(**hole)
        elif isinstance(getattr(cfg,name),C.Array):getattr(cfg,name)[:]=value
        else:setattr(cfg,name,value)
    if list(cfg.n)!=record['resolution']:raise ValueError('checkpoint/configuration grid mismatch')
    n,_,nphi=record['resolution']
    with np.load(ROOT/f"validation/raw/{record['case']}_{n}_{nphi}.npz") as saved:
        unknowns=saved['unknowns'].copy()
    return cfg,unknowns
