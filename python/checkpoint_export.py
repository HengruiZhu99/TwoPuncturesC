"""Portable, explicit-field HiSpID sampler checkpoint (text format version1).

This is a data interchange file, not a raw ctypes structure. Seventeen digit
values preserve IEEE double unknowns. Library SHA and acceptance describe the
source evidence; the consumer must explicitly name that SHA in its input.
"""
from pathlib import Path
import hashlib
import numpy as np
from hispid import Config

PARAMETERIZATION='W_plus_Aminus1_V'

def write_checkpoint(path,config,unknowns,library_sha256,acceptance):
    if len(library_sha256)!=64 or any(c not in '0123456789abcdef' for c in library_sha256):
        raise ValueError('invalid library SHA256')
    if acceptance not in ('analytic_seed','preliminary','strong','diagnostic'):
        raise ValueError('invalid acceptance label')
    values=np.asarray(unknowns,dtype=float).ravel()
    if values.size!=4*np.prod(list(config.n)) or not np.all(np.isfinite(values)):
        raise ValueError('invalid checkpoint unknowns')
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='ascii',newline='\n') as out:
        out.write('HISPID_CHECKPOINT 1\nparameterization '+PARAMETERIZATION+'\n')
        out.write('library_sha256 '+library_sha256+'\nacceptance '+acceptance+'\n')
        for name,_ in Config._fields_:
            value=getattr(config,name)
            if name=='hole':
                for index,hole in enumerate(value):
                    numbers=[hole.mass,*hole.center,*hole.spin,*hole.velocity]
                    out.write('hole'+str(index)+' '+' '.join(format(x,'.17g') for x in numbers)+'\n')
            else:
                numbers=list(value) if hasattr(value,'_length_') else [value]
                out.write(name+' '+' '.join(format(x,'.17g') if isinstance(x,float) else str(x) for x in numbers)+'\n')
        out.write('unknowns '+str(values.size)+'\n')
        for start in range(0,len(values),4):
            out.write(' '.join(format(x,'.17g') for x in values[start:start+4])+'\n')
        out.write('END\n')
    return dict(path=str(path.resolve()),file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                source_library_sha256=library_sha256,acceptance=acceptance)
