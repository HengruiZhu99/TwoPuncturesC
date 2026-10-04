"""Analytic controls for the independent Decimal four-metric oracle."""
from decimal import Decimal, localcontext
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'validation'))
from boosted_schwarzschild_oracle import BoostedSchwarzschild, decimal_input, sample_float64, compare_retained


def relative(a, b):
    return abs(a-b)/(1+abs(b))


class BoostedSchwarzschildOracleTests(unittest.TestCase):
    def test_binary64_inputs_are_exact(self):
        self.assertEqual(decimal_input(.1), Decimal.from_float(.1))
        self.assertNotEqual(decimal_input(.1), Decimal('.1'))
        for value in (float('nan'), float('inf'), Decimal('NaN'), Decimal('Infinity')):
            with self.assertRaises(ValueError):decimal_input(value)

    def test_unboosted_static_and_flat_controls(self):
        for m, v in ((1, (0, 0, 0)), (0, ('.2', '-.3', '.1'))):
            oracle = BoostedSchwarzschild(m, v)
            point = ('1.1', '.2', '-.3')
            f = oracle.fields(point)
            self.assertTrue(all(x == 0 for row in f['Kij'] for x in row))
            self.assertTrue(all(x == 0 for x in f['shift']))
            if not m:
                self.assertEqual(f['gamma'], [[Decimal(i == j) for j in range(3)] for i in range(3)])
                self.assertLess(abs(f['lapse']-1), Decimal('1e-65'))
            c = oracle.constraints(point)
            self.assertLess(abs(c['H']), Decimal('1e-65'))
            self.assertEqual(c['Mnorm'], 0)

    def test_regular_K_agrees_with_four_christoffel_definition(self):
        oracle = BoostedSchwarzschild(1, ('.2', '-.3', '.1'))
        for point in (('1.1', '.2', '.3'), ('.08', '.03', '-.02')):
            regular = oracle.fields(point)['Kij']; direct = oracle.direct_adm_K(point)
            with localcontext() as ctx:
                ctx.prec = 80
                self.assertLess(max(relative(a, b) for ra, rb in zip(regular, direct) for a, b in zip(ra, rb)), Decimal('1e-65'))

    def test_signed_lapse_exact_throat_and_vacuum_constraints(self):
        oracle = BoostedSchwarzschild(1, ('.6', '0', '0'))  # G=1.25 exactly.
        self.assertLess(oracle.fields(('.2', '0', '0'))['lapse'], 0)
        self.assertEqual(oracle.fields(('.4', '0', '0'))['lapse'], 0)
        self.assertGreater(oracle.fields(('.8', '0', '0'))['lapse'], 0)
        with self.assertRaises(ValueError):oracle.direct_adm_K(('.4', '0', '0'))
        for point in (('.4', '0', '0'), ('.2', '.01', '-.03'), ('.8', '.02', '.03')):
            c = oracle.constraints(point)
            self.assertLess(abs(c['H']), Decimal('1e-62'))
            self.assertLess(c['Mnorm'], Decimal('1e-39'))

    def test_generic_metric_derivatives_and_reversed_boost(self):
        oracle = BoostedSchwarzschild(1, ('.2', '-.3', '.1'))
        reverse = BoostedSchwarzschild(1, ('-.2', '.3', '-.1'))
        point = [Decimal('1.1'), Decimal('.2'), Decimal('.3')]
        f = oracle.fields(point); b = reverse.fields(point)
        with localcontext() as ctx:
            ctx.prec = 80
            self.assertEqual(f['gamma'], b['gamma'])
            self.assertLess(max(abs(a+b) for ra, rb in zip(f['Kij'], b['Kij']) for a, b in zip(ra, rb)), Decimal('1e-65'))
            h = Decimal('1e-20')
            for d in range(3):
                left = point.copy(); right = point.copy();left[d] -= h;right[d] += h
                a = oracle.fields(left)['gamma'];b = oracle.fields(right)['gamma']
                self.assertLess(max(relative((b[i][j]-a[i][j])/(2*h), f['dgamma'][d][i][j])
                    for i in range(3) for j in range(3)), Decimal('1e-37'))

    def test_gamma10_near_horizon_and_precision_independence(self):
        v = (0.99498743710662, 0., 0.)
        point = (.075, .002, -.003)
        a = BoostedSchwarzschild(1., v, precision=70)
        b = BoostedSchwarzschild(1., v, precision=90)
        fa, fb = a.fields(point), b.fields(point)
        with localcontext() as ctx:
            ctx.prec = 90
            for key in ('gamma', 'Kij'):
                self.assertLess(max(relative(x, y) for ra, rb in zip(fa[key], fb[key]) for x, y in zip(ra, rb)), Decimal('1e-60'))
        for oracle in (a, b):
            c = oracle.constraints(point)
            self.assertLess(abs(c['H']), Decimal('1e-55'))
            self.assertLess(c['Mnorm'], Decimal('1e-34'))

    def test_reject_invalid_physical_inputs(self):
        with self.assertRaises(ValueError):BoostedSchwarzschild(-1)
        with self.assertRaises(ValueError):BoostedSchwarzschild(1, (1, 0, 0))
        with self.assertRaises(ValueError):BoostedSchwarzschild(1).fields((0, 0, 0))
        with self.assertRaises(ValueError):BoostedSchwarzschild(1).fields((1, 2, 3, 4))

    def test_retained_derivative_layout_and_unchanged_field_bound(self):
        import numpy as np
        oracle = BoostedSchwarzschild(1., (.6, 0., 0.))
        xyz = np.array([[.8, .02, .03]], dtype=np.float64)
        values = sample_float64(oracle, xyz)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'witness.npz'
            for shape in ((1, 27), (1, 3, 3, 3)):
                values['dgamma'] = values['dgamma'].reshape(shape)
                np.savez(path, xyz=xyz, **values)
                result = compare_retained(oracle, path, [0])
                self.assertTrue(result['passed'])
                self.assertEqual(result['field_scaled_linf_limit'], 1e-12)
                self.assertFalse(result['physical_acceptance'])
            values['Kij'][0, 0] += 1e-8
            np.savez(path, xyz=xyz, **values)
            self.assertFalse(compare_retained(oracle, path, [0])['passed'])


if __name__ == '__main__':
    unittest.main()
