"""Verify explicit native library images, including dyld install-name reuse."""
import ctypes as C
import hashlib
from pathlib import Path
import sys

_loaded_builds={}

def verify_image(lib,symbol,expected):
    expected=Path(expected).resolve(strict=True)
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
        if actual!=expected:
            raise ValueError(f'native loader reused {actual}; compare distinct builds in separate processes')
    digest=hashlib.sha256(expected.read_bytes()).hexdigest()
    previous=_loaded_builds.get(address)
    if previous is not None and previous!=(expected,digest):
        raise ValueError('loaded native library differs from on-disk build; start a fresh process')
    _loaded_builds[address]=(expected,digest)
    return digest

def validate_krylov(krylov,linear_rtol):
    if krylov not in (None,'gmres','bicgstab'):raise ValueError('krylov must be gmres or bicgstab')
    if linear_rtol is not None:
        import math
        if not math.isfinite(linear_rtol) or not 0<linear_rtol<1:
            raise ValueError('linear_rtol must be finite and between 0 and 1; None selects native forcing')
