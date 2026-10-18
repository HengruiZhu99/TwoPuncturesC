"""Checked budget rejection before allocating a large BY modal factor bank."""
import argparse,ctypes as C,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();lib=C.CDLL(str(a.library.resolve()))
f=lib.TP_modal_create_analytic_with_budget;f.argtypes=[C.c_int]*3+[C.POINTER(C.c_double),C.c_int];f.restype=C.c_void_p
old=lib.TP_modal_create_analytic;old.argtypes=[C.c_int]*3+[C.POINTER(C.c_double)];old.restype=C.c_void_p
u=np.zeros(64*128*16);ptr=u.ctypes.data_as(C.POINTER(C.c_double));cases=[]
for na,nb,np_,budget in [(64,128,16,16),(12,18,8,-1),(12,18,8,15),(12,18,8,65537),(2**31-1,2**31-1,256,65536)]:
 assert not f(na,nb,np_,ptr,budget)
 cases.append([na,nb,np_,budget])
# The old public entry point must still reject a >4GiB bank. This call
# cannot read U: the checked capacity rejection precedes all seed access.
assert not old(224,448,16,ptr)
a.output.write_text(json.dumps(dict(passed=True,rejected_cases=cases,legacy_four_gib_guard_preserved=True),indent=2)+'\n');print('checked budget/overflow guards pass')
