"""Direct axisymmetric single-hole grids, independent of the two-focus map.

CPU tensor products are economical for these small two-dimensional problems.
The BY class solves the original vacuum Hamiltonian equation. The Kerr class
interpolates exact native seed geometry; it is not an evolution or a solver
for perturbed Kerr free data. Both use r=L*A/(1-A), 0<A<1, and mu=cos(theta).
"""
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, gmres
from numpy.polynomial.legendre import leggauss, legvander


class SingleHoleGrid:
    def __init__(self, nr, ntheta, scale):
        if not isinstance(nr, int) or not isinstance(ntheta, int) or not 4 <= nr <= 512 or not 4 <= ntheta <= 256:
            raise ValueError('grid requires 4<=nr<=512 and 4<=ntheta<=256')
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError('positive finite radial scale required')
        self.nr, self.ntheta, self.scale = nr, ntheta, float(scale)
        t = np.pi*(np.arange(nr)+.5)/nr
        self.A = (1-np.cos(t))/2
        self.mu, w = leggauss(ntheta)
        self.weights = (-1.)**np.arange(nr)*np.sin(t)
        dx = self.A[:, None]-self.A[None, :]
        np.fill_diagonal(dx, 1.)
        D = self.weights[None, :]/self.weights[:, None]/dx
        np.fill_diagonal(D, 0.)
        np.fill_diagonal(D, -D.sum(axis=1))
        self.D = D
        self.D2 = D@D
        self.L = legvander(self.mu, ntheta-1)
        self.to_modes = self.L.T*w[None, :]*(2*np.arange(ntheta)+1)[:, None]/2

    def coordinates(self):
        r = self.scale*self.A/(1-self.A)
        r, mu = np.broadcast_arrays(r[:, None], self.mu[None, :])
        return np.stack((r*np.sqrt(1-mu*mu), np.zeros_like(r), r*mu), axis=-1)

    def evaluate(self, values, radius, mu):
        radius, mu = np.broadcast_arrays(np.asarray(radius, float), np.asarray(mu, float))
        if np.any(~np.isfinite(radius)) or np.any(radius <= 0) or np.any(~np.isfinite(mu)) or np.any(abs(mu) > 1+1e-13):
            raise ValueError('finite positive radii and |mu|<=1 required')
        return self.evaluate_compactified(values, radius/(radius+self.scale), mu)

    def evaluate_compactified(self, values, A, mu):
        values = np.asarray(values, dtype=float)
        if values.shape[:2] != (self.nr, self.ntheta) or not np.all(np.isfinite(values)):
            raise ValueError('finite values on this grid are required')
        A, mu = np.broadcast_arrays(np.asarray(A, float), np.asarray(mu, float))
        if np.any(~np.isfinite(A)) or np.any((A < 0) | (A > 1)) or np.any(~np.isfinite(mu)) or np.any(abs(mu) > 1+1e-13):
            raise ValueError('compactified coordinates require 0<=A<=1 and |mu|<=1')
        shape = A.shape
        A = A.reshape(-1)
        dx = A[:, None]-self.A[None, :]
        exact = dx == 0
        safe = np.where(exact, 1, dx)
        W = self.weights[None, :]/safe
        W /= W.sum(axis=1)[:, None]
        for p in np.flatnonzero(exact.any(axis=1)):
            W[p] = exact[p]
        angular = legvander(np.clip(mu.reshape(-1), -1, 1), self.ntheta-1)
        coefficients = np.einsum('lm,im...->il...', self.to_modes, values)
        result = np.einsum('pi,pl,il...->p...', W, angular, coefficients)
        return result.reshape(shape+values.shape[2:])


class BowenYorkSingleHole:
    """Aligned spin/linear momentum along z; no fictitious second puncture.

    S and P are BY parameters, not measured horizon spin and velocity.
    V=u/(1-A) remains bounded at both compactified ends. The preconditioner
    inverts the exact separated flat operator in Legendre modes.
    """
    def __init__(self, mass=1., spin=0., momentum=0., nr=32, ntheta=24):
        if not np.isfinite(mass) or mass <= 0 or not np.isfinite(spin) or not np.isfinite(momentum):
            raise ValueError('finite BY charges and positive mass required')
        self.mass, self.spin, self.momentum = float(mass), float(spin), float(momentum)
        self.grid = g = SingleHoleGrid(nr, ntheta, mass/2)
        A = g.A
        self.radial = (A*A*(1-A)**2)[:, None]*g.D2+(2*A*(1-A)*(1-2*A))[:, None]*g.D
        self.radial -= np.diag(2*A*(1-A))
        self.ell = np.arange(ntheta)*(np.arange(ntheta)+1)
        self.factors = [lu_factor(self.radial-np.diag(np.full(nr, l))) for l in self.ell]
        self.V = np.zeros((nr, ntheta))

    def linear(self, V):
        g = self.grid
        return self.radial@V-(V@g.to_modes.T*self.ell)@g.L.T

    def inverse(self, rhs):
        g = self.grid
        modes = rhs@g.to_modes.T
        return np.stack([lu_solve(f, modes[:, l]) for l, f in enumerate(self.factors)], axis=1)@g.L.T

    def source(self, V):
        g = self.grid
        A, mu = g.A[:, None], g.mu[None, :]
        denominator = 1+A*(1-A)*V
        if np.any(denominator <= 0):
            return np.full_like(V, np.nan), np.full_like(V, np.nan)
        numerator = ((9/4)*self.spin**2/g.scale**4*A**3*(1-A)**3*(1-mu*mu)
                     +(9/16)*self.momentum**2/g.scale**2*A**5*(1-A)*(1+2*mu*mu))
        source = numerator/denominator**7
        return source, -7*A*(1-A)*source/denominator

    def solve(self, tolerance=2e-13, max_newton=20):
        if not np.isfinite(tolerance) or tolerance <= 0 or not isinstance(max_newton, int) or max_newton < 1:
            raise ValueError('positive tolerance and Newton budget required')
        shape, size = self.V.shape, self.V.size
        steps = 0
        for iteration in range(max_newton+1):
            source, slope = self.source(self.V)
            residual = self.linear(self.V)+source
            norm = float(np.max(abs(residual)))
            if np.isfinite(norm) and norm <= tolerance:
                return dict(converged=True, newton_iterations=iteration, krylov_iterations=steps, scaled_linf=norm)
            if iteration == max_newton or not np.isfinite(norm):
                break
            operator = LinearOperator((size, size), matvec=lambda x:(self.linear(x.reshape(shape))+slope*x.reshape(shape)).ravel())
            preconditioner = LinearOperator((size, size), matvec=lambda x:self.inverse(x.reshape(shape)).ravel())
            count = [0]
            def callback(_):
                count[0] += 1
            change, status = gmres(operator, -residual.ravel(), M=preconditioner, rtol=1e-5, atol=0., restart=40, maxiter=8, callback=callback, callback_type='pr_norm')
            steps += count[0]
            if status:
                break
            for backtrack in range(12):
                trial = self.V+2.**(-backtrack)*change.reshape(shape)
                new_source, _ = self.source(trial)
                trial_norm = np.max(abs(self.linear(trial)+new_source))
                if np.isfinite(trial_norm) and trial_norm < norm:
                    self.V = trial
                    break
            else:
                break
        return dict(converged=False, newton_iterations=iteration, krylov_iterations=steps, scaled_linf=norm)

    def regularized_correction(self, radius, mu):
        return self.grid.evaluate(self.V, radius, mu)

    def sample(self, xyz):
        x = np.asarray(xyz, float).reshape(-1, 3)
        r = np.linalg.norm(x, axis=1)
        if np.any(r <= 0):
            raise ValueError('puncture is excluded')
        n = x/r[:, None]
        u = self.grid.scale/(r+self.grid.scale)*self.regularized_correction(r, n[:, 2])
        psi = 1+self.mass/(2*r)+u
        if np.any(~np.isfinite(psi)) or np.any(psi <= 0):
            raise ValueError('nonpositive interpolated conformal factor; refine the grid')
        P, S = np.array([0., 0., self.momentum]), np.array([0., 0., self.spin])
        outer = lambda a, b:np.einsum('ni,nj->nij', a, b)
        p = np.broadcast_to(P, n.shape)
        A = 1.5/r[:, None, None]**2*(outer(p, n)+outer(n, p)
            -(np.eye(3)-outer(n, n))*(n@P)[:, None, None])
        cross = np.cross(S, n)
        A += 3/r[:, None, None]**3*(outer(cross, n)+outer(n, cross))
        return dict(gamma=psi[:, None, None]**4*np.broadcast_to(np.eye(3), A.shape),
                    Kij=A/psi[:, None, None]**2, psi=psi, attenuation=np.ones(len(x)))


class KerrSingleHole:
    """Spectral representation of Kerr with spin and optional boost along z.

    Native analytic geometry is evaluated only at grid construction. Sampling
    subsequently interpolates regularized spherical tensor components. Thus
    refinement tests measure representation error, not repeated exact calls.
    """
    def __init__(self, backend, mass=1., chi=.99, nr=32, ntheta=24, seed_family='liu', velocity=0.):
        from hispid import Hole
        if not np.isfinite(mass) or mass <= 0 or not np.isfinite(chi) or abs(chi) >= 1:
            raise ValueError('positive mass and subextremal spin required')
        if seed_family not in ('qi', 'liu'):
            raise ValueError('supported Kerr charts are qi and liu')
        if not np.isfinite(velocity) or abs(velocity) >= 1:
            raise ValueError('subluminal boost required')
        self.mass, self.chi, self.seed_family = float(mass), float(chi), seed_family
        self.velocity = float(velocity)
        self.grid = g = SingleHoleGrid(nr, ntheta, mass/2)
        self.horizon_radius = mass*(1+np.sqrt(1-chi*chi))/4 if seed_family == 'liu' else mass*np.sqrt(1-chi*chi)/2
        xyz = g.coordinates()
        r = np.linalg.norm(xyz, axis=-1)
        mu = xyz[..., 2]/r
        sint = np.sqrt(1-mu*mu)
        er = xyz/r[..., None]
        et = np.stack((mu, np.zeros_like(mu), -sint), axis=-1)
        ep = np.zeros_like(er)
        ep[..., 1] = 1
        frame = np.stack((er, et, ep), axis=-1)
        exact = backend.seed(Hole(mass, spin=(0, 0, chi*mass*mass), velocity=(0, 0, velocity)), xyz.reshape(-1, 3), seed_family=seed_family)
        gamma, K = [frame.swapaxes(-1, -2)@exact[k].reshape(nr, ntheta, 3, 3)@frame for k in ('gamma', 'Kij')]
        A = g.A[:, None]
        weight_g = A**4/(1-A)
        self.k_inner_power, self.k_outer_power = (4, 2) if velocity else (1, 3)
        weight_k = A**self.k_inner_power/(1-A)**self.k_outer_power
        def regularize(tensor, identity):
            return np.stack((tensor[..., 0, 0]-identity, tensor[..., 1, 1]-identity,
                (tensor[..., 2, 2]-tensor[..., 1, 1])/(1-mu*mu),
                tensor[..., 0, 1]/sint, tensor[..., 0, 2]/sint,
                tensor[..., 1, 2]/(1-mu*mu)), axis=-1)
        self.values = np.concatenate((weight_g[..., None]*regularize(gamma, 1),
                                     weight_k[..., None]*regularize(K, 0)), axis=-1)

    def sample(self, xyz):
        x = np.asarray(xyz, float).reshape(-1, 3)
        r = np.linalg.norm(x, axis=1)
        if np.any(r <= 0):
            raise ValueError('puncture is excluded')
        er = x/r[:, None]
        mu, sint = er[:, 2], np.hypot(er[:, 0], er[:, 1])
        A = r/(r+self.grid.scale)
        v = self.grid.evaluate(self.values, r, mu)
        w = np.stack((-er[:, 1], er[:, 0], np.zeros(len(x))), axis=-1)
        et = np.zeros_like(er)
        nonaxis = sint > 0
        et[nonaxis, :2] = er[nonaxis, :2]*mu[nonaxis, None]/sint[nonaxis, None]
        et[:, 2] = -sint
        outer = lambda a, b:np.einsum('ni,nj->nij', a, b)
        def reconstruct(values, identity):
            rr, tt, difference, rt, rp, tp = values.T
            result = (identity+tt)[:, None, None]*np.eye(3)+(rr-tt)[:, None, None]*outer(er, er)
            result += difference[:, None, None]*outer(w, w)
            result += (rt*sint)[:, None, None]*(outer(er, et)+outer(et, er))
            result += rp[:, None, None]*(outer(er, w)+outer(w, er))
            result += (tp*sint)[:, None, None]*(outer(et, w)+outer(w, et))
            return result
        gamma = reconstruct(((1-A)/A**4)[:, None]*v[:, :6], 1.)
        K = reconstruct(((1-A)**self.k_outer_power/A**self.k_inner_power)[:, None]*v[:, 6:], 0.)
        if not np.all(np.isfinite(gamma)) or not np.all(np.isfinite(K)) or np.any(np.linalg.eigvalsh(gamma) <= 0):
            raise ValueError('nonpositive or nonfinite interpolated Kerr geometry; refine the grid')
        return dict(gamma=gamma, Kij=K, psi=np.linalg.det(gamma)**(1/12), attenuation=np.ones(len(x)))

    def horizon(self, order=64):
        if self.velocity:
            raise ValueError('this spherical horizon diagnostic requires a rest hole')
        mu, weight = leggauss(order)
        sint = np.sqrt(1-mu*mu)
        r = self.horizon_radius
        er = np.stack((sint, np.zeros_like(mu), mu), axis=-1)
        et = np.stack((mu, np.zeros_like(mu), -sint), axis=-1)
        def components(radius):
            q = self.sample(radius*er)
            contract = lambda v:np.einsum('ni,nij,nj->n', v, q['gamma'], v)
            return q, contract(er), contract(et), q['gamma'][:, 1, 1]
        q, gr, gt, gp = components(r)
        area = 2*np.pi*r*r*np.sum(weight*np.sqrt(gt*gp))
        krp = np.einsum('ni,ni->n', er, q['Kij'][:, :, 1])
        J = r**3/4*np.sum(weight*krp*sint*np.sqrt(gt*gp/gr))
        Mirr = np.sqrt(area/(16*np.pi))
        mass = np.sqrt(Mirr*Mirr+J*J/(4*Mirr*Mirr))
        h = r*2e-4
        logs = []
        for shift in (-2, -1, 1, 2):
            _, _, t, p = components(r+shift*h)
            logs.append(.5*np.log(t*p))
        derivative = np.einsum('a,an->n', np.array([1, -8, 8, -1])/12, logs)/h
        expansion = (2/r+derivative)/np.sqrt(gr)
        return dict(area=float(area), angular_momentum=float(J), irreducible_mass=float(Mirr),
                    mass=float(mass), chi=float(J/mass**2), expansion_rms=float(np.sqrt(np.mean(expansion**2))))
