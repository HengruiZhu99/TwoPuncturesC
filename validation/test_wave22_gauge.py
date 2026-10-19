"""Regression: reconstruct the high-spin packet with its actual Kerr gauge."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'validation')]
from hispid import Backend,Wave22
from kerr_wave_case import config
from wave22_study import multipoles


def run(b):
    chi=.9999;threshold=chi/(1+np.sqrt(1-chi*chi))
    w=Wave22(amplitude=.001,omega=.99*threshold,cycles=4.,r0=12.,cutoff_min=2.,cutoff_max=4.)
    radius=w.r0+np.pi*w.cycles/w.omega;y=radius-w.r0
    W=2*np.tanh(y)
    h=np.sqrt(12)*w.amplitude*W/radius*np.exp(-1j*(w.omega*y+np.pi/4))
    expected=dict(metric=h,dt_metric=-1j*w.omega*h,dr_metric=(-1/radius-1j*w.omega)*h)
    c=config(b,chi=chi);c.hole[1].center[:]=(0,0,-32)
    with b.create_sampler(c) as seed,b.create_sampler(c,wave=w) as trial:
        baseline=multipoles(seed,radius,4,chi);actual=multipoles(trial,radius,4,chi)
        errors={key:abs(complex(*actual[key])-complex(*baseline[key])-value) for key,value in expected.items()}
        assert max(errors.values())<1e-12,errors
        wrong=multipoles(trial,radius,4,.9)
        wrong_error=abs(complex(*wrong['dt_metric'])-complex(*baseline['dt_metric'])-expected['dt_metric'])
        assert wrong_error>100*max(errors.values()),(wrong_error,errors)
    return dict(passed=True,chi=chi,seed_superradiant_threshold=threshold,wave22=w.as_dict(),
                multipole_radius=radius,absolute_multipole_errors=errors,wrong_spin_gauge_dt_error=wrong_error,
                library_sha256=b.library_sha256())

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path)
    a=p.parse_args();result=run(Backend(str(a.library.resolve())));print(json.dumps(result,indent=2))
    if a.output:a.output.write_text(json.dumps(result,indent=2)+'\n')
