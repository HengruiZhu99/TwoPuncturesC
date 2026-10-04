"""Independent high-precision isotropic Schwarzschild boost oracle.

No native code, Kerr formulas, graph-slice construction, or Jet implementation
is used. Units are G=c=1, signature is (-,+,+,+), and
K_ij = -(partial_t gamma_ij - D_i beta_j - D_j beta_i)/(2 alpha).
The rest coordinates are T=G(t-v.x), X=Lx-Gvt, where
L=I+(G-1)vv^T/v^2. Thus the black hole moves with lab velocity +v.

In rest coordinates p=(1+m/(2r))^4 and a=(r-m/2)/(r+m/2). Writing
b_mu=(G,-Gv) and q=p-a^2 gives g_mu_nu=p eta_mu_nu+q b_mu b_nu.
The exterior lapse is positive. Its signed continuation across the isotropic
throat is used below; using a positive lapse on both sheets flips inner K.
The four-metric is coordinate-degenerate exactly at r=m/2, but the regular
spatial metric and K formula extend there by continuity.

Ricci uses analytic Cartesian metric derivatives. Momentum uses independent
fourth-order Decimal differences of K with a tiny configurable step. These
are supplemental isolated-seed checks, never binary acceptance or a replacement
for the existing failed binary64 Cartesian finite-difference gates.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path

ZERO, ONE, TWO = Decimal(0), Decimal(1), Decimal(2)
FIELD_TOLERANCE = 1e-12  # Existing portable physical-field/gradient bound.


def decimal_input(value):
    """Preserve binary64 input bits; strings instead denote exact decimals."""
    if isinstance(value, Decimal):
        if not value.is_finite():
            raise ValueError('finite oracle inputs required')
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError('finite oracle inputs required')
        return Decimal.from_float(value)
    result = Decimal(value)
    if not result.is_finite():
        raise ValueError('finite oracle inputs required')
    return result


def _matrix(n):
    return [[ZERO for _ in range(n)] for _ in range(n)]


class BoostedSchwarzschild:
    """Fixed physical parameters; arithmetic uses local Decimal contexts."""

    def __init__(self, mass, velocity=(0, 0, 0), center=(0, 0, 0), precision=80):
        self.mass = decimal_input(mass)
        self.velocity = tuple(decimal_input(x) for x in velocity)
        self.center = tuple(decimal_input(x) for x in center)
        if self.mass < 0 or len(self.velocity) != 3 or len(self.center) != 3:
            raise ValueError('nonnegative mass and three-vector parameters required')
        if not isinstance(precision, int) or precision < 40:
            raise ValueError('at least 40 decimal digits required')
        self.precision = precision
        with localcontext() as ctx:
            ctx.prec = precision
            self.v2 = sum(x*x for x in self.velocity)
            if self.v2 >= 1:
                raise ValueError('subluminal velocity required')
            self.lorentz = ONE/(ONE-self.v2).sqrt()
            self.spatial_boost = tuple(tuple(
                Decimal(i == j)+(self.lorentz-ONE)*self.velocity[i]*self.velocity[j]/self.v2
                if self.v2 else Decimal(i == j) for j in range(3)) for i in range(3))

    def _state(self, point, time=0):
        point = tuple(point)
        if len(point) != 3:
            raise ValueError('three Cartesian coordinates required')
        x = tuple(decimal_input(q)-c for q, c in zip(point, self.center))
        t = decimal_input(time)
        v, G, L = self.velocity, self.lorentz, self.spatial_boost
        y = [sum(L[i][j]*x[j] for j in range(3))-G*v[i]*t for i in range(3)]
        r = sum(q*q for q in y).sqrt()
        if r == 0:
            if self.mass:
                raise ValueError('isotropic puncture excluded')
            r = ONE  # Flat metric is constant, including at the coordinate origin.
        J = [[-G*v[i], *L[i]] for i in range(3)]
        dr = [sum(y[k]*J[k][mu] for k in range(3))/r for mu in range(4)]
        ddr = [[(sum(J[k][mu]*J[k][nu] for k in range(3))-dr[mu]*dr[nu])/r
                for nu in range(4)] for mu in range(4)]
        z = self.mass/(TWO*r)
        p = (ONE+z)**4
        a = (ONE-z)/(ONE+z)
        pp = -4*z*(ONE+z)**3/r
        ppp = (12*z*z*(ONE+z)**2+8*z*(ONE+z)**3)/(r*r)
        ap = TWO*z/(r*(ONE+z)**2)
        app = -4*z/(r*r*(ONE+z)**3)
        q, qp, qpp = p-a*a, pp-TWO*a*ap, ppp-TWO*(ap*ap+a*app)
        b = [G, *(-G*w for w in v)]
        eta = [Decimal(-1), ONE, ONE, ONE]
        four = [[p*eta[i]*Decimal(i == j)+q*b[i]*b[j] for j in range(4)] for i in range(4)]
        prime = [[pp*eta[i]*Decimal(i == j)+qp*b[i]*b[j] for j in range(4)] for i in range(4)]
        second = [[ppp*eta[i]*Decimal(i == j)+qpp*b[i]*b[j] for j in range(4)] for i in range(4)]
        g = [row[1:] for row in four[1:]]
        dg = [[[prime[i+1][j+1]*dr[k+1] for j in range(3)] for i in range(3)] for k in range(3)]
        ddg = [[[[second[i+1][j+1]*dr[k+1]*dr[l+1]+prime[i+1][j+1]*ddr[k+1][l+1]
                  for j in range(3)] for i in range(3)] for l in range(3)] for k in range(3)]
        D = p-a*a*self.v2
        if D <= 0:
            raise ValueError('nonspacelike slice')
        inv = [[(Decimal(i == j)-q*v[i]*v[j]/D)/p for j in range(3)] for i in range(3)]
        C = p.sqrt()/(G*D.sqrt())
        alpha = a*C
        beta = [-q*w/D for w in v]
        # Contracting K_ij=-alpha Gamma^0_ij cancels every inverse a:
        # sum_sigma g^{0sigma} partial_sigma r=-partial_t r/p.
        K = [[G*G*C*(ap-a*pp/(TWO*p))*(v[i]*dr[j+1]+v[j]*dr[i+1])
              -a*C*dr[0]*prime[i+1][j+1]/(TWO*p) for j in range(3)] for i in range(3)]
        return dict(gamma=g, Kij=K, dgamma=dg, ddgamma=ddg, inverse=inv,
                    lapse=alpha, shift=beta, rest_radius=r, rest_lapse=a,
                    four_metric=four, radial_derivative=prime, radial_gradient=dr,
                    p=p, q=q, b=b, eta=eta)

    def fields(self, point, time=0):
        """Decimal physical gamma, K, Cartesian dgamma, signed lapse and shift."""
        with localcontext() as ctx:
            ctx.prec = self.precision
            s = self._state(point, time)
            return {k: s[k] for k in ('gamma', 'Kij', 'dgamma', 'lapse', 'shift', 'rest_radius', 'rest_lapse')}

    def direct_adm_K(self, point, time=0):
        """Uncancelled four-Christoffel definition, independent of regular K."""
        with localcontext() as ctx:
            ctx.prec = self.precision
            s = self._state(point, time)
            a, p, q, b, eta = s['rest_lapse'], s['p'], s['q'], s['b'], s['eta']
            if not a:
                raise ValueError('four-metric chart degenerates at exact throat')
            upper_b = [eta[i]*b[i] for i in range(4)]
            inv_time = [eta[0]*Decimal(mu == 0)/p-q*upper_b[0]*upper_b[mu]/(p*a*a) for mu in range(4)]
            prime, dr = s['radial_derivative'], s['radial_gradient']
            return [[-s['lapse']/TWO*sum(inv_time[mu]*(prime[mu][j+1]*dr[i+1]
                +prime[mu][i+1]*dr[j+1]-prime[i+1][j+1]*dr[mu]) for mu in range(4))
                for j in range(3)] for i in range(3)]

    def constraints(self, point, step=None):
        """High-precision H and M using analytic metric derivatives and FD K.

        M_cov[i]=D_j K^j_i-partial_i trace(K); M is its raised version.
        The K derivative step is separate from the existing binary64 verifier.
        """
        with localcontext() as ctx:
            ctx.prec = self.precision
            x = tuple(decimal_input(q) for q in point)
            h = decimal_input(step) if step is not None else max(self.mass, ONE)*Decimal('1e-12')
            if h <= 0:
                raise ValueError('positive momentum derivative step required')
            s = self._state(x)
            g, inv, dg, ddg, K = (s[k] for k in ('gamma', 'inverse', 'dgamma', 'ddgamma', 'Kij'))
            dinv = [[[-sum(inv[i][a]*dg[d][a][b]*inv[b][j] for a in range(3) for b in range(3))
                      for j in range(3)] for i in range(3)] for d in range(3)]
            conn = [[[sum(inv[k][l]*(dg[i][l][j]+dg[j][l][i]-dg[l][i][j]) for l in range(3))/TWO
                      for j in range(3)] for i in range(3)] for k in range(3)]
            dc = [[[[sum(dinv[d][k][l]*(dg[i][l][j]+dg[j][l][i]-dg[l][i][j])
                  +inv[k][l]*(ddg[d][i][l][j]+ddg[d][j][l][i]-ddg[d][l][i][j])
                  for l in range(3))/TWO for j in range(3)] for i in range(3)] for k in range(3)] for d in range(3)]
            Ric = [[sum(dc[k][k][i][j]-dc[j][k][i][k]
                    +sum(conn[k][i][j]*conn[l][k][l]-conn[l][i][k]*conn[k][j][l] for l in range(3))
                    for k in range(3)) for j in range(3)] for i in range(3)]
            scalar = sum(inv[i][j]*Ric[i][j] for i in range(3) for j in range(3))
            mixed = [[sum(inv[i][a]*K[a][j] for a in range(3)) for j in range(3)] for i in range(3)]
            trace = sum(mixed[i][i] for i in range(3))
            square = sum(mixed[i][j]*mixed[j][i] for i in range(3) for j in range(3))
            dk = []
            for d in range(3):
                values = []
                for offset in (-2, -1, 1, 2):
                    y = list(x); y[d] += offset*h
                    values.append(self._state(y)['Kij'])
                dk.append([[sum(Decimal(w)*value[i][j] for w, value in zip((1, -8, 8, -1), values))/(12*h)
                            for j in range(3)] for i in range(3)])
            dt = [sum(dinv[d][a][b]*K[a][b]+inv[a][b]*dk[d][a][b]
                      for a in range(3) for b in range(3)) for d in range(3)]
            momentum_low = [sum(sum(dinv[j][j][a]*K[a][i]+inv[j][a]*dk[j][a][i] for a in range(3))
                +sum(conn[j][j][k]*mixed[k][i]-conn[k][j][i]*mixed[j][k] for k in range(3))
                for j in range(3))-dt[i] for i in range(3)]
            momentum = [sum(inv[i][j]*momentum_low[j] for j in range(3)) for i in range(3)]
            norm2 = sum(momentum_low[i]*inv[i][j]*momentum_low[j] for i in range(3) for j in range(3))
            return dict(H=scalar+trace*trace-square, M=momentum, M_cov=momentum_low,
                        Mnorm=max(ZERO, norm2).sqrt(), R=scalar, trace=trace, K2=square,
                        precision=self.precision, momentum_derivative_step=h)


def sample_float64(oracle, xyz):
    """Round exact oracle fields once for the unchanged physical.py FD oracle."""
    import numpy as np
    points = np.asarray(xyz, dtype=np.float64).reshape(-1, 3)
    if not len(points):
        raise ValueError('at least one Cartesian point required')
    values = [oracle.fields(tuple(float(x) for x in point)) for point in points]
    return {key: np.asarray([value[key] for value in values], dtype=np.float64).reshape(len(points), -1)
            for key in ('gamma', 'Kij', 'dgamma')} | {'attenuation': np.ones(len(points))}


def compare_retained(oracle, witness, indices):
    """Selected immutable physical-field checks only; provenance belongs to caller."""
    import numpy as np
    path = Path(witness).resolve(strict=True)
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    with np.load(path, allow_pickle=False) as saved:
        xyz = saved['xyz']
        if xyz.dtype != np.float64 or xyz.ndim != 2 or xyz.shape[1] != 3:
            raise ValueError('binary64 Cartesian witness coordinates required')
        if not indices or len(set(indices)) != len(indices) or any(i < 0 or i >= len(xyz) for i in indices):
            raise ValueError('distinct in-range witness indices required')
        rows = []
        for index in indices:
            point = tuple(float(x) for x in xyz[index])
            exact = oracle.fields(point)
            differences = {}
            for key, width in (('gamma', 9), ('Kij', 9), ('dgamma', 27)):
                if saved[key].dtype != np.float64 or saved[key].shape != (len(xyz), width):
                    raise ValueError('binary64 witness shape differs: '+key)
                native = saved[key][index]
                if not np.isfinite(native).all():
                    raise ValueError('nonfinite retained fields')
                flat = np.asarray(exact[key], dtype=object).ravel()
                with localcontext() as ctx:
                    ctx.prec = oracle.precision
                    errors = [abs(decimal_input(float(a))-b)/(ONE+abs(b)) for a, b in zip(native, flat)]
                differences[key] = float(max(errors))
            rows.append(dict(index=index, xyz=list(point), differences=differences,
                             passed=all(q <= FIELD_TOLERANCE for q in differences.values())))
    if hashlib.sha256(path.read_bytes()).hexdigest() != before:
        raise ValueError('witness changed while reading')
    return dict(witness=str(path), witness_sha256=before, rows=rows,
                field_scaled_linf_limit=FIELD_TOLERANCE, passed=all(r['passed'] for r in rows),
                purpose='supplemental_selected_isolated_fields_only', physical_acceptance=False,
                supersedes_failed_seed_controls=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mass', type=float, default=1.)
    p.add_argument('--velocity', nargs=3, type=float, default=[math.sqrt(.99), 0., 0.])
    p.add_argument('--center', nargs=3, type=float, default=[0., 0., 0.])
    p.add_argument('--point', nargs=3, type=float, default=[.075, .002, -.003])
    p.add_argument('--precision', type=int, default=80)
    p.add_argument('--momentum-step', default='1e-12')
    p.add_argument('--witness', help='retained producer NPZ; no native library is loaded')
    p.add_argument('--indices', default='0', help='explicit comma-separated retained point indices')
    a = p.parse_args()
    oracle = BoostedSchwarzschild(a.mass, a.velocity, a.center, a.precision)
    result = dict(parameters=dict(mass=a.mass, velocity=a.velocity, center=a.center,
        binary64_inputs_exact=True, lorentz=str(oracle.lorentz)),
        point=a.point, fields=oracle.fields(a.point), constraints=oracle.constraints(a.point, a.momentum_step),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), physical_acceptance=False)
    if a.witness:
        result['retained_comparison'] = compare_retained(oracle, a.witness, list(map(int, a.indices.split(','))))
    print(json.dumps(result, default=lambda x: str(x) if isinstance(x, Decimal) else x, indent=2))
    return int(a.witness is not None and not result['retained_comparison']['passed'])


if __name__ == '__main__':
    raise SystemExit(main())
