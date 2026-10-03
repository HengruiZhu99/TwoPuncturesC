"""Explicit execution-space selection and fenced native work accounting."""
import ctypes as C
import json

class Statistics(C.Structure):
    _fields_=[(name,C.c_double) for name in ('operator_seconds','precondition_seconds','setup_seconds','linear_seconds','transfer_seconds','nonlinear_seconds')]+[(name,C.c_ulonglong) for name in ('resident_bytes','peak_workspace_bytes','host_to_device_bytes','device_to_host_bytes','kokkos_host_current_bytes','kokkos_host_peak_bytes','kokkos_device_current_bytes','kokkos_device_peak_bytes','nonlinear_calls','workspace_creations')]+[('memory_tracking_available',C.c_int)]

class BYStatistics(C.Structure):
    _fields_=[(name,C.c_int) for name in ('workspace_creations','device_residual_calls','original_confirmation_calls','polishing_steps')]+[(name,C.c_double) for name in ('device_candidate_linf','original_confirmation_linf')]

def by_statistics(lib):
    if not hasattr(lib,'TP_solver_get_execution_statistics'):return None
    lib.TP_solver_get_execution_statistics.argtypes=[C.POINTER(BYStatistics)];lib.TP_solver_get_execution_statistics.restype=None
    s=BYStatistics();lib.TP_solver_get_execution_statistics(C.byref(s));return {key:getattr(s,key) for key,_ in s._fields_}

def select(lib,execution='reference',threads=0):
    if execution not in ('reference','kokkos'):raise ValueError('execution must be reference or kokkos')
    if execution=='reference':return 0
    if not hasattr(lib,'Puncture_execution_initialize'):raise ValueError('library lacks Kokkos execution API')
    lib.Puncture_execution_initialize.argtypes=[C.c_int];lib.Puncture_execution_initialize.restype=C.c_int
    if lib.Puncture_execution_initialize(threads):raise ValueError('Kokkos unavailable or initialization conflicts with the current runtime')
    return 1

def name(lib):
    if not hasattr(lib,'Puncture_execution_name'):return 'reference'
    lib.Puncture_execution_name.restype=C.c_char_p
    return lib.Puncture_execution_name().decode('ascii')

def concurrency(lib):
    if not hasattr(lib,'Puncture_execution_concurrency'):return 1
    lib.Puncture_execution_concurrency.restype=C.c_int
    return lib.Puncture_execution_concurrency()

def device_description(lib):
    if not hasattr(lib,'Puncture_execution_device_description'):return None
    lib.Puncture_execution_device_description.argtypes=[C.c_char_p,C.c_int];lib.Puncture_execution_device_description.restype=C.c_int
    buffer=C.create_string_buffer(1024)
    if lib.Puncture_execution_device_description(buffer,len(buffer)):return None
    return json.loads(buffer.value)

def statistics(lib):
    if not hasattr(lib,'Puncture_execution_statistics'):return None
    lib.Puncture_execution_statistics.argtypes=[C.POINTER(Statistics)];lib.Puncture_execution_statistics.restype=C.c_int
    s=Statistics()
    if lib.Puncture_execution_statistics(C.byref(s)):raise ValueError('execution statistics unavailable')
    return {k:getattr(s,k) for k,_ in s._fields_}

def reset_statistics(lib):
    if hasattr(lib,'Puncture_execution_reset_statistics'):
        lib.Puncture_execution_reset_statistics.argtypes=[];lib.Puncture_execution_reset_statistics.restype=None
        lib.Puncture_execution_reset_statistics()
