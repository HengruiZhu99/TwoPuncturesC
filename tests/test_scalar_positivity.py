"""Closed-form reconstruction controls for the limited scalar diagnosis."""
import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'validation'))
from inspect_scalar_positivity import scalar_correction,seed_scalar


class ScalarReconstruction(unittest.TestCase):
    def test_constant_and_quadratic_cartesian_modes(self):
        x=np.array([[1.3,.2,.4],[-.8,.7,-.3],[3.5,-1.,.6]])
        b=2.;shape=(12,16,8);grid=np.zeros((*shape[::-1],4))
        total=sum(np.linalg.norm(x-np.array([sign*b,0,0]),axis=1) for sign in (-1,1))
        t=(total-2*b)/(total+2*b);phi=np.arctan2(x[:,2],x[:,1]);rho=np.hypot(x[:,1],x[:,2])
        for mode in (0,2):
            grid[:]=0;grid[mode,:,:,0]=1
            norm=np.sqrt((1 if mode==0 else 2)/shape[2])
            expected=-2*(1-t)*(rho*(1-t)/(2*b))**mode*norm*np.cos(mode*phi)
            for lam,kappa in ((.2,2.),(.03,4.),(1.,.5)):
                actual=scalar_correction(x,grid.ravel(),shape,dict(radial_stretch=lam,angular_stretch=kappa),b)
                np.testing.assert_allclose(actual,expected,rtol=1e-13,atol=1e-14)

    def test_zero_spin_QI_reduces_to_isotropic_scalar(self):
        holes=[dict(mass=.5,spin=[0.,0.,0.],center=[sign*2.,0.,0.]) for sign in (-1,1)]
        x=np.array([[1.3,.2,.4],[-.8,.7,-.3],[3.5,-1.,.6]])
        expected=1+sum(h['mass']/(2*np.linalg.norm(x-np.asarray(h['center']),axis=1)) for h in holes)
        np.testing.assert_allclose(seed_scalar(x,dict(hole=holes)),expected,rtol=1e-14,atol=1e-14)


if __name__=='__main__':unittest.main()
