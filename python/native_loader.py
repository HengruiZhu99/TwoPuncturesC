"""Verify explicit native library images, including dyld install-name reuse."""
import ctypes as C
import hashlib
from pathlib import Path
import sys

_loaded_builds={}

def image_path(lib,symbol):
    address=C.cast(getattr(lib,symbol),C.c_void_p).value
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
        return actual
    raise ValueError('native image verification requires dladdr')

def verify_image(lib,symbol,expected):
    expected=Path(expected).resolve(strict=True)
    address=C.cast(getattr(lib,symbol),C.c_void_p).value
    actual=image_path(lib,symbol)
    if actual!=expected:
        raise ValueError(f'native loader reused {actual}; compare distinct builds in separate processes')
    digest=hashlib.sha256(expected.read_bytes()).hexdigest()
    previous=_loaded_builds.get(address)
    if previous is not None and previous!=(expected,digest):
        raise ValueError('loaded native library differs from on-disk build; start a fresh process')
    _loaded_builds[address]=(expected,digest)
    return digest

def verify_puncture_dependencies(lib,primary):
    primary=Path(primary).resolve(strict=True);images={}
    for symbol in ('AB_To_XR','C_To_c','PK_solve','Puncture_execution_initialize'):
        if not hasattr(lib,symbol):continue
        actual=image_path(lib,symbol)
        if actual!=primary and (actual.parent!=primary.parent or actual.name not in ('libTwoPunctures.so','libTwoPunctures.dylib')):
            raise ValueError(f'unexpected puncture dependency for {symbol}: {actual}')
        images[str(actual)]=verify_image(lib,symbol,actual)
    return images

def loaded_kokkos_images():
    """Record the shared runtime actually mapped, rather than a build glob."""
    paths=set()
    if sys.platform.startswith('linux'):
        for line in Path('/proc/self/maps').read_text().splitlines():
            parts=line.split(maxsplit=5)
            if len(parts)==6 and parts[5].startswith('/') and 'libkokkos' in Path(parts[5]).name:
                paths.add(Path(parts[5]).resolve(strict=True))
    elif sys.platform=='darwin':
        dyld=C.CDLL(None);dyld._dyld_image_count.restype=C.c_uint
        dyld._dyld_get_image_name.argtypes=[C.c_uint];dyld._dyld_get_image_name.restype=C.c_char_p
        for i in range(dyld._dyld_image_count()):
            name=dyld._dyld_get_image_name(i)
            if name and 'libkokkos' in Path(name.decode()).name:paths.add(Path(name.decode()).resolve(strict=True))
    else:raise ValueError('runtime image verification requires Linux or macOS')
    return {str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(paths)}

def validate_krylov(krylov,linear_rtol,allow_lgmres=False):
    if krylov not in ((None,'gmres','bicgstab','lgmres') if allow_lgmres else (None,'gmres','bicgstab')):raise ValueError('krylov must be gmres or bicgstab')
    if linear_rtol is not None:
        import math
        if not math.isfinite(linear_rtol) or not 0<linear_rtol<1:
            raise ValueError('linear_rtol must be finite and between 0 and 1; None selects native forcing')
