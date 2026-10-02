"""Choose an equation system independently of the linear Krylov backend.

Backend('hispid', path) / Backend('bowen_york', path) expose config(), create().
Both solutions support solve(krylov='gmres'|'bicgstab', linear_rtol=...).
Each system keeps its own physical configuration and preconditioner options.
"""
def Backend(system,library):
    if system=='hispid':
        from hispid import Backend as Implementation
    elif system=='bowen_york':
        from bowen_york import Backend as Implementation
    else:raise ValueError('system must be hispid or bowen_york')
    return Implementation(library)
