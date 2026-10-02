"""Public system/backend selector and option/lifetime guards on small data."""
import os
from pathlib import Path
import unittest
import numpy as np
from punctures import Backend
from hispid import SolveOptions
import ctypes as C

class Selectors(unittest.TestCase):
    def test_invalid_system(self):
        with self.assertRaises(ValueError):Backend('other','/unused')
    def test_hi_options_and_backends(self):
        b=Backend('hispid',os.environ['HISPID_LIBRARY']);c=b.config()
        c.n[:]=(12,18,8);c.far_radius=0;c.max_krylov=2000;c.tolerance=1e-13
        c.hole[0].spin[2]=.03;c.hole[1].velocity[1]=-.02
        results=[]
        for method in ('gmres','bicgstab'):
            with b.create(c) as s:
                for invalid in (0,-1,1,float('nan')):
                    with self.assertRaises(ValueError):s.solve(krylov=method,linear_rtol=invalid)
                with self.assertRaises(ValueError):s.solve(krylov='bad')
                o=SolveOptions(0,0,.001)
                self.assertEqual(b.lib.HiSpID_solve_with_options(s.context,C.byref(o)),-1)
                r=s.solve(krylov=method,linear_rtol=.001)
                self.assertEqual(r['status'],0,r)
                self.assertEqual(s.resolved_options['krylov'],method)
                self.assertTrue(all(row[2]<=.001 for row in s.linear_history()))
                self.assertTrue(np.isfinite(s.unknowns()).all())
                results.append(s.sample([[0,2,1],[.2,.3,.4],[5,3,2]]))
        # Interface control checks physical fields. The separate matrix report
        # retains the failed raw-P/coefficient equivalence gate.
        for key in ('gamma','Kij','psi','correction'):
            self.assertLess(np.max(np.abs(results[0][key]-results[1][key])/(1+np.abs(results[0][key]))),1e-10)
    def test_by_options_lifetime_and_backends(self):
        hi=Path(os.environ['HISPID_LIBRARY']);path=hi.parent.parent/'lib/libTwoPunctures.so'
        if not path.exists():self.skipTest('BY shared library not built')
        b=Backend('bowen_york',str(path));c=b.config()
        c['integer'].update(npoints_A=8,npoints_B=12,npoints_phi=8)
        c['real'].update(par_S_plus3=.03,par_P_minus2=-.008)
        results=[]
        for method in ('gmres','bicgstab'):
            with b.create(c) as s:
                with self.assertRaises(ValueError):b.create(c)
                with self.assertRaises(ValueError):s.sample([[0,2,1]])
                with self.assertRaises(ValueError):s.solve(krylov=method,linear_rtol=0)
                r=s.solve(krylov=method,linear_rtol=.001,preconditioner='modal',max_krylov=2000)
                self.assertEqual(r['status'],0,r)
                self.assertLessEqual(r['residual_linf'],1e-12)
                self.assertEqual(s.work_statistics()['linear_failures'],0)
                results.append(s.unknowns())
                self.assertTrue(np.isfinite(s.sample([[0,2,1]])['gamma']).all())
                with self.assertRaises(ValueError):s.solve()
            with self.assertRaises(ValueError):s.sample([[0,2,1]])
        self.assertLess(np.max(np.abs(results[0]-results[1])/(1+np.abs(results[0]))),1e-10)

        c['integer']['Newton_maxit']=0
        with b.create(c) as s:
            r=s.solve(krylov='gmres',preconditioner='modal')
            self.assertFalse(r['converged']);self.assertEqual(r['status'],1)
        c['integer']['Newton_maxit']=24;c['integer']['TP_krylov_solver']=9
        with b.create(c) as s:
            c_bad=s.config['integer'];c_bad['TP_preconditioner']=9
            # Reserved solver keys are intentionally resolved by solve().
            s.config['real']['Newton_tol']='bad'
            with self.assertRaises(ValueError):s.solve()
            s.config['real']['Newton_tol']=1e-12
            self.assertEqual(s.solve(preconditioner='modal')['status'],0)
