"""Portable, explicit-field HiSpID sampler checkpoint (text format version1).

This is a data interchange file, not a raw ctypes structure. Seventeen digit
values preserve IEEE double unknowns. Library SHA and acceptance describe the
source evidence; the consumer must explicitly name that SHA in its input.
"""
from pathlib import Path
import hashlib
import numpy as np
import re
from hispid import Config,seed_family_code
import ctypes as C

PARAMETERIZATION='modal_P_C2prolate_mapped_v2'

def _validate_parameterization(parameterization):
    custom_map=re.fullmatch(r'modal_P_(?:C2prolate_map_v3|C4prolate_map_v4|C2tauC4_map_v5)_r([0-9eE+.-]+)_k([0-9eE+.-]+)',parameterization)
    if parameterization not in ('W_plus_Aminus1_V','W_plus_Aminus1_V_C2axis','modal_P_C2prolate_mapped_v2') and custom_map is None:
        raise ValueError('unsupported checkpoint parameterization')
    if custom_map and not (.001<=float(custom_map[1])<=1 and .1<=float(custom_map[2])<=6):
        raise ValueError('unsupported collocation maps')

def read_checkpoint(path):
    """Read a current typed checkpoint without loading any native image.

    This only decodes data. The caller must verify the producer SHA, continuous
    basis/maps, and any explicit sampler migration before using its values.
    """
    path=Path(path).resolve(strict=True)
    initial_sha=hashlib.sha256(path.read_bytes()).hexdigest()
    with path.open(encoding='ascii') as stream:
        magic=stream.readline().split()
        if magic not in (['HISPID_CHECKPOINT','1'],['HISPID_CHECKPOINT','2']):
            raise ValueError('unsupported checkpoint version')
        version=int(magic[1])
        header={}
        for line in stream:
            words=line.split()
            if not words or words[0] in header:raise ValueError('empty or duplicate checkpoint field')
            if words[0]=='unknowns':
                if len(words)!=2:raise ValueError('invalid unknown count')
                count=int(words[1]);break
            header[words[0]]=words[1:]
        else:raise ValueError('missing checkpoint unknowns')
        expected={'parameterization','library_sha256','acceptance','hole0','hole1'} | {name for name,_ in Config._fields_ if name!='hole'}
        if version==2:expected.add('seed_family')
        family='qi'
        if version==2:
            if len(header.get('seed_family',[]))!=1:raise ValueError('missing seed family')
            family=header['seed_family'][0];seed_family_code(family)
        if set(header)!=expected:raise ValueError('checkpoint configuration inventory differs')
        if any(len(header[k])!=1 for k in ('parameterization','library_sha256','acceptance')):
            raise ValueError('invalid checkpoint metadata')
        sha=header['library_sha256'][0]
        if re.fullmatch('[0-9a-f]{64}',sha) is None:raise ValueError('invalid producer SHA256')
        if header['acceptance'][0] not in ('analytic_seed','preliminary','strong','diagnostic'):
            raise ValueError('invalid checkpoint acceptance')
        _validate_parameterization(header['parameterization'][0])
        config=Config();config.seed_family=family
        for name,kind in Config._fields_:
            if name=='hole':
                for h in range(2):
                    numbers=list(map(float,header['hole'+str(h)]))
                    if len(numbers)!=10 or not np.isfinite(numbers).all():raise ValueError('invalid hole data')
                    config.hole[h].mass=numbers[0]
                    for field,start in (('center',1),('spin',4),('velocity',7)):
                        getattr(config.hole[h],field)[:]=numbers[start:start+3]
                continue
            value=getattr(config,name)
            subtype=kind._type_ if isinstance(value,C.Array) else kind
            numbers=[int(x) if subtype is C.c_int else float(x) for x in header[name]]
            if len(numbers)!=(len(value) if isinstance(value,C.Array) else 1) or not np.isfinite(numbers).all():
                raise ValueError('invalid checkpoint field: '+name)
            if subtype is C.c_int and any(not -2**31<=x<2**31 for x in numbers):
                raise ValueError('checkpoint integer outside native range: '+name)
            if isinstance(value,C.Array):value[:]=numbers
            else:setattr(config,name,numbers[0])
        if min(config.n)<4 or count!=4*int(config.n[0])*int(config.n[1])*int(config.n[2]) or not 16<=config.memory_limit_mib<=65536:
            raise ValueError('invalid checkpoint grid/budget/count')
        if (count//4)*128>config.memory_limit_mib*1024**2:
            raise ValueError('checkpoint exceeds minimum sampler budget')
        values=[]
        for line in stream:
            if line.split()==['END']:break
            values.extend(map(float,line.split()))
            if len(values)>count:raise ValueError('extra checkpoint unknowns')
        else:raise ValueError('missing checkpoint terminator')
        if stream.read().strip():raise ValueError('trailing checkpoint data')
    values=np.asarray(values,dtype=float)
    if values.size!=count or not np.isfinite(values).all():raise ValueError('invalid checkpoint unknowns')
    if hashlib.sha256(path.read_bytes()).hexdigest()!=initial_sha:raise ValueError('checkpoint changed while decoding')
    return config,values,dict(path=str(path),file_sha256=initial_sha,
        source_library_sha256=sha,acceptance=header['acceptance'][0],parameterization=header['parameterization'][0],seed_family=family,format_version=version)

def write_checkpoint(path,config,unknowns,library_sha256,acceptance,parameterization=PARAMETERIZATION,seed_family=None):
    _validate_parameterization(parameterization)
    configured=getattr(config,'seed_family','qi')
    if seed_family is not None and hasattr(config,'seed_family') and seed_family!=configured:
        raise ValueError('seed-family override conflicts with configuration metadata')
    family=configured if seed_family is None else seed_family
    seed_family_code(family);version=1 if family=='qi' else 2
    if len(library_sha256)!=64 or any(c not in '0123456789abcdef' for c in library_sha256):
        raise ValueError('invalid library SHA256')
    if acceptance not in ('analytic_seed','preliminary','strong','diagnostic'):
        raise ValueError('invalid acceptance label')
    values=np.asarray(unknowns,dtype=float).ravel()
    if values.size!=4*np.prod(list(config.n)) or not np.all(np.isfinite(values)):
        raise ValueError('invalid checkpoint unknowns')
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='ascii',newline='\n') as out:
        out.write('HISPID_CHECKPOINT '+str(version)+'\nparameterization '+parameterization+'\n')
        if version==2:out.write('seed_family '+family+'\n')
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
                source_library_sha256=library_sha256,acceptance=acceptance,seed_family=family,format_version=version)
