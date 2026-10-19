"""Independent analytic controls for the optional direct spherical path."""
import sys
import unittest
import os
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'python'))
from single_hole import SingleHoleGrid, BowenYorkSingleHole


class SingleHoleTests(unittest.TestCase):
    def test_tensor_interpolation_and_endpoints(self):
        g = SingleHoleGrid(16, 12, .5)
        A, mu = g.A[:, None], g.mu[None, :]
        values = A**3*(3*mu*mu-1)/2
        a = np.array([0., .13, .72, 1.])
        m = np.array([-1., -.2, .6, 1.])
        self.assertLess(np.max(abs(g.evaluate_compactified(values, a, m)-a**3*(3*m*m-1)/2)), 2e-13)
        self.assertLess(np.max(abs(g.D@g.A**3-3*g.A**2)), 2e-12)
        self.assertLess(np.max(abs(g.D2@g.A**3-6*g.A)), 2e-10)
        # Collocation nodes must be reproduced, not divided by zero.
        np.testing.assert_allclose(g.evaluate_compactified(values, A, mu), values, rtol=0, atol=2e-13)

    def test_brill_lindquist_and_small_spin(self):
        s = BowenYorkSingleHole()
        self.assertTrue(s.solve()['converged'])
        x = np.array([[1., .2, .3], [.1, -.2, .07]])
        np.testing.assert_allclose(s.sample(x)['psi'], 1+.5/np.linalg.norm(x, axis=1), rtol=0, atol=1e-14)
        s = BowenYorkSingleHole(spin=.001)
        self.assertTrue(s.solve(tolerance=1e-15)['converged'])
        A, mu = np.array([.2, .4, .7, .9]), np.array([.1, .3, .5, .8])
        u = (1-A)*s.grid.evaluate_compactified(s.V, A, mu)/1e-6
        # Leading small-spin solution of the flat Hamiltonian equation.
        # The quadrupole is negative for the standard P2=(3mu²-1)/2.
        analytic = .4*(-2*A**5+6*A**4-5*A**3+1)-.8*(1-A)**3*A*A*(3*mu*mu-1)/2
        self.assertLess(np.max(abs(u-analytic)), 5e-7)

    def test_invalid_configuration(self):
        for args in ((3, 12, .5), (12, 3, .5), (12, 12, 0.)):
            with self.assertRaises(ValueError):
                SingleHoleGrid(*args)
        for kw in ({'mass':0.}, {'spin':np.nan}, {'momentum':np.inf}):
            with self.assertRaises(ValueError):
                BowenYorkSingleHole(**kw)

    @unittest.skipUnless(os.environ.get('HISPID_TEST_LIBRARY'), 'explicit native image required')
    def test_native_tensor_and_modal_replay(self):
        from hispid import Backend, Hole
        sys.path.insert(0, str(ROOT/'validation'))
        from binary_reference_norm import modal_U
        b = Backend(os.environ['HISPID_TEST_LIBRARY'])
        x = np.array([[.4, .2, .3], [2., .4, -.5], [0., 0., 1.], [0., 0., -2.]])
        k = b.single_hole(mass=1.7, chi=.4, nr=64, ntheta=24)
        exact = b.seed(Hole(1.7, spin=(0, 0, .4*1.7**2)), x, seed_family='liu')
        for name in ('gamma', 'Kij'):
            expected = exact[name].reshape(-1, 3, 3)
            self.assertLess(np.max(abs(k.sample(x)[name]-expected)/(1+abs(expected))), 2e-9)
        k = b.single_hole(chi=0., velocity=.2, nr=32, ntheta=24, seed_family='qi')
        exact = b.seed(Hole(1., velocity=(0, 0, .2)), x)
        for name in ('gamma', 'Kij'):
            expected = exact[name].reshape(-1, 3, 3)
            self.assertLess(np.max(abs(k.sample(x)[name]-expected)/(1+abs(expected))), 2e-9)
        c = b.config()
        c.n[:] = [12, 24, 8]
        c.far_radius = 0
        coords = np.array([[.2, .3, .4], [.7, -.4, 1.2]])
        A, B, phi = coords.T
        a2 = A*A
        bhalf = (c.hole[0].center[0]-c.hole[1].center[0])/2
        x = np.stack((bhalf*(a2+1)/(a2-1)*2*B/(1+B*B),
            bhalf*2*A/(1-a2)*(1-B*B)/(1+B*B)*np.cos(phi),
            bhalf*2*A/(1-a2)*(1-B*B)/(1+B*B)*np.sin(phi)), axis=-1)
        with b.create_sampler(c) as s:
            values = np.random.default_rng(123).normal(size=s.size)*1e-4
            s.set_unknowns(values)
            observed = s.sample(x)['correction']/(A-1)[:, None]
            np.testing.assert_allclose(modal_U(values, list(c.n), coords), observed, rtol=1e-10, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
