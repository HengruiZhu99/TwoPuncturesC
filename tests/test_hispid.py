"""Fast scientifically meaningful controls. Run with explicit HISPID_LIBRARY."""
import os
import unittest
import numpy as np
from hispid import Backend,Hole
from physical import constraints,norms,charges,extrapolate
from prolong import prolong

class HiSpIDTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.b=Backend(os.environ['HISPID_LIBRARY'])
        cls.x=np.array([[1.1,.2,.3],[.7,-.9,.4],[2.7,.5,-1.2]])

    def test_schwarzschild_and_conformal_choices(self):
        h=Hole(1);v=self.b.seed(h,self.x);p=1+.5/np.linalg.norm(self.x,axis=1)
        np.testing.assert_allclose(v['psi'],p,rtol=1e-14)
        np.testing.assert_allclose(v['gamma'].reshape(-1,3,3),p[:,None,None]**4*np.eye(3),rtol=1e-14,atol=1e-14)
        np.testing.assert_array_equal(v['Kij'],0)
        h=Hole(1,spin=(.2,.3,.4),velocity=(.2,-.1,.15))
        a=self.b.seed(h,self.x,0);b=self.b.seed(h,self.x,1)
        np.testing.assert_allclose(a['gamma'],b['gamma'],rtol=1e-14)
        np.testing.assert_allclose(a['Kij'],b['Kij'],rtol=1e-14)
        np.testing.assert_allclose(np.linalg.det(b['conformal_metric'].reshape(-1,3,3)),1,rtol=1e-14)

    def test_grid_prolongation_polynomial_and_nyquist(self):
        def values(shape):
            a,b,p=shape
            A=-np.cos(np.pi*(np.arange(a)+.5)/a)[None,None,:]
            B=-np.cos(np.pi*(np.arange(b)+.5)/b)[None,:,None]
            phi=2*np.pi*np.arange(p)[:,None,None]/p
            v=.1*A*A*B+.02*B**3*np.cos(2*phi)+.03*A*B*np.sin(phi)+.04*np.cos(4*phi)
            return (v[:,:,:,None]*np.array([1.,.4,-.2,.3])).ravel()
        old=(8,7,8);new=(13,12,12)
        np.testing.assert_allclose(prolong(values(old),old,new),values(new),atol=1e-14,rtol=1e-14)

    def test_physical_seed_step_convergence(self):
        for spin,velocity in [((0,0,0),(0,0,0)),((0,0,.6),(0,0,0)),((.2,.3,.4),(.2,-.1,.15))]:
            sample=lambda x:self.b.seed(Hole(1,spin=spin,velocity=velocity),x)
            results=[norms(constraints(sample,self.x,step)) for step in (.03,.015,.0075)]
            self.assertLess(results[-1]['H_rms'],1e-8)
            self.assertGreater(results[0]['H_rms']/results[1]['H_rms'],10)
            self.assertGreater(results[1]['H_rms']/results[2]['H_rms'],10)
            self.assertLess(results[-1]['M_rms'],1e-8)

    def test_seed_covariance(self):
        axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis);angle=.64
        A=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
        Q=np.eye(3)+np.sin(angle)*A+(1-np.cos(angle))*(A@A);offset=np.array([.37,-.21,.13])
        spin=np.array([.2,.3,.4]);vel=np.array([.2,-.1,.15])
        a=self.b.seed(Hole(1,spin=spin,velocity=vel),self.x)
        b=self.b.seed(Hole(1,center=offset,spin=Q@spin,velocity=Q@vel),self.x@Q.T+offset)
        for key in ('gamma','Kij','conformal_metric','Atilde'):
            expected=np.einsum('ik,nkl,jl->nij',Q,a[key].reshape(-1,3,3),Q).reshape(-1,9)
            np.testing.assert_allclose(b[key],expected,rtol=1e-12,atol=1e-13)

    def test_boosted_schwarzschild_signed_mean_curvature(self):
        # Paper II.B.2, including inner sheet r<m/2; fixes lapse orientation.
        m=1.;v=.4;G=1/np.sqrt(1-v*v);h=Hole(m,velocity=(0,v,0))
        x=np.array([[.9,.3,.2],[.08,.09,.07],[.44,.1,.08]])
        r=np.sqrt(x[:,0]**2+(G*x[:,1])**2+x[:,2]**2)
        B=np.sqrt((m+2*r)**6-16*(m-2*r)**2*r**4*v*v)
        expected=32*G*m*v*((m+2*r)**7-32*(m-2*r)**2*(m-r)*r**4*v*v)*r**2*x[:,1]/((m+2*r)**3*B**3)
        out=self.b.seed(h,x,0)
        np.testing.assert_allclose(out['mean_curvature'],expected,rtol=1e-11,atol=1e-13)

    def test_analytic_jvp_and_native_solve(self):
        c=self.b.config();c.n[:]=[8,8,6];c.hole[0].spin[:]=[0,0,.05];c.hole[1].spin[:]=[.02,0,0]
        c.hole[0].velocity[:]=[.03,.05,0];c.hole[1].velocity[:]=[-.01,-.03,.02]
        c.inner_max[:]=[0,0];c.far_radius=0
        with self.b.create(c) as s:
            rng=np.random.default_rng(14108607);base=rng.normal(0,1e-6,s.size);d=rng.normal(0,1e-3,s.size)
            epsilon=1e-4;fd=(s.residual(base+epsilon*d)-s.residual(base-epsilon*d))/(2*epsilon)
            exact=s.jvp(base,d)
            self.assertLess(np.linalg.norm(exact-fd)/np.linalg.norm(exact),1e-8)
            diag=s.solve();self.assertEqual(diag['status'],0,diag)
            self.assertLess(max(diag['scaled_linf']),c.tolerance)
            self.assertGreater(np.max(np.abs(s.sample(self.x)['correction'][:,1:])),1e-7)

    def test_brill_lindquist_limit_and_contexts(self):
        c=self.b.config();c.n[:]=[6,6,4];c.conformal_choice=0;c.far_radius=0;c.inner_max[:]=[0,0]
        with self.b.create(c) as s,self.b.create(c) as t:
            self.assertEqual(s.solve()['status'],0)
            self.assertEqual(t.solve()['status'],0)
            out=s.sample(self.x);p=np.ones(len(self.x))
            for hole in c.hole:p+=hole.mass/(2*np.linalg.norm(self.x-np.array(hole.center),axis=1))
            np.testing.assert_allclose(out['psi'],p,rtol=1e-14)
            self.assertLess(np.max(np.abs(out['correction'])),1e-14)
            q=s.charges(100,ntheta=8,nphi=16)
            self.assertAlmostEqual(q[0],(1+1/200)**3,delta=.002)
        with self.assertRaises(ValueError):s.sample(self.x)

    def test_nyquist_vector_gradient_and_distant_axis(self):
        # The cosine Nyquist mode has zero d/dphi at collocation points,
        # but its analytic derivative must survive off-grid reconstruction.
        c=self.b.config();c.n[:]=[6,6,8];c.conformal_choice=0
        c.far_radius=0;c.inner_max[:]=[0,0];c.inner_flatten=0
        c.hole[0]=Hole(.5,(1,0,0));c.hole[1]=Hole(.5,(-1,0,0))
        m=4;x=np.array([[.3,.7,.4],[1.7,-.8,.2],[-.4,.3,-.9]])
        with self.b.create(c) as s:
            u=np.zeros((8,6,6,4))
            for k in range(8):u[k,:,:,1]=np.cos(m*2*np.pi*k/8)
            s.set_unknowns(u.ravel())
            rp=np.linalg.norm(x-np.array([1.,0,0]),axis=1)
            rm=np.linalg.norm(x+np.array([1.,0,0]),axis=1)
            X=np.arccosh((rp+rm)/2);A=2*np.tanh(X/2)-1
            dX=((x-np.array([1.,0,0]))/rp[:,None]+(x+np.array([1.,0,0]))/rm[:,None])/(2*np.sinh(X)[:,None])
            dA=dX/np.cosh(X/2)[:,None]**2
            phi=np.arctan2(x[:,2],x[:,1]);rho2=x[:,1]**2+x[:,2]**2
            dphi=np.c_[np.zeros(len(x)),-x[:,2]/rho2,x[:,1]/rho2]
            grad=dA*np.cos(m*phi)[:,None]-(A-1)[:,None]*m*np.sin(m*phi)[:,None]*dphi
            L=np.zeros((len(x),3,3));L[:,0,:]+=grad;L[:,:,0]+=grad
            L-=2/3*grad[:,0,None,None]*np.eye(3)
            np.testing.assert_allclose(s.sample(x)['Atilde'].reshape(-1,3,3),L,rtol=1e-12,atol=1e-12)
            s.set_unknowns(np.zeros(s.size))
            axis=s.sample([[1e5,0,0],[-1e5,0,0],[0,0,0]])
            self.assertTrue(np.all(np.isfinite(axis['gamma'])))

    def test_analytic_far_split_preserves_brill_lindquist(self):
        # An exact solution with a far filter requires u=sum((1-F)m/2r).
        # The fixed analytic part should carry it even on a tiny grid.
        c=self.b.config();c.n[:]=[6,6,4];c.conformal_choice=0
        c.inner_max[:]=[0,0];c.inner_flatten=0;c.far_radius=40
        x=np.r_[self.x,[[30,11,-8],[80,-20,17],[400,50,-60]]]
        with self.b.create(c) as s:
            d=s.solve();self.assertEqual(d['status'],0,d)
            self.assertLess(d['newton_iterations'],2)
            out=s.sample(x);psi=np.ones(len(x));W=np.zeros(len(x))
            for hole in c.hole:
                r=np.linalg.norm(x-np.array(hole.center),axis=1)
                psi+=hole.mass/(2*r);W+=-np.expm1(-(r/c.far_radius)**4)*hole.mass/(2*r)
            np.testing.assert_allclose(out['psi'],psi,rtol=2e-14)
            np.testing.assert_allclose(out['correction'][:,0],W,atol=2e-14)

    def test_far_split_jvp(self):
        c=self.b.config();c.n[:]=[6,6,6];c.far_radius=8
        c.hole[0].spin[:]=[.04,.01,.03];c.hole[0].velocity[:]=[.04,.02,-.03]
        with self.b.create(c) as s:
            rng=np.random.default_rng(17);base=rng.normal(0,1e-5,s.size);direction=rng.normal(0,.01,s.size)
            eps=1e-4
            fd=(s.residual(base+eps*direction)-s.residual(base-eps*direction))/(2*eps)
            exact=s.jvp(base,direction)
            self.assertLess(np.linalg.norm(fd-exact)/np.linalg.norm(exact),1e-8)

    def test_far_split_retains_unattenuated_laplacian(self):
        c=self.b.config();c.n[:]=[6,6,4];c.conformal_choice=0;c.far_radius=2
        c.hole[0]=Hole(1,(3,0,0));c.hole[1]=Hole(0,(-3,0,0))
        c.inner_min[:]=[.1,0];c.inner_max[:]=[.4,0]
        with self.b.create(c) as s:
            v=s.equation_samples();r=np.linalg.norm(v['xyz']-np.array(c.hole[0].center),axis=1)
            F=np.exp(-(r/2)**4);lapW=.5*F*(12*r/2**4-16*r**5/2**8)
            expected=-8*(1-v['attenuation'])*lapW/v['psi']**5
            self.assertGreater(np.sum(v['attenuation']<1),0)
            self.assertGreater(np.max(np.abs(expected)),1e-5)
            np.testing.assert_allclose(v['physical_equivalent'][:,0],expected,atol=1e-12)

    def test_domain_and_invalid_input(self):
        with self.assertRaises(ValueError):self.b.seed(Hole(1),[[0,0,0]])
        throat=self.b.seed(Hole(1),[[.5,0,0]])
        np.testing.assert_array_equal(throat['Kij'],0)
        np.testing.assert_allclose(throat['gamma'],np.eye(3).reshape(1,9)*16)
        c=self.b.config();c.hole[0].velocity[:]=[1,0,0]
        with self.assertRaises(ValueError):self.b.create(c)
        c=self.b.config();c.n[:]=[256,256,256]
        with self.assertRaises(ValueError):self.b.create(c)
        c=self.b.config();c.n[:]=[4,4,4]
        with self.b.create(c) as s:
            with self.assertRaises(ValueError):s.charges(10,center=[0,0])
        with self.assertRaises(ValueError):self.b.operators(c,[0,0],np.zeros(40))

    def test_regular_boosted_kerr_throat(self):
        spin=np.array([.2,.3,.4]);velocity=np.array([.4,-.2,.15]);G=1/np.sqrt(1-velocity@velocity)
        rh=.5*np.sqrt(1-spin@spin);direction=np.array([1.,2.,-3.]);direction/=np.linalg.norm(direction)
        rest=rh*direction
        lab=rest+(1/G-1)*velocity*(velocity@rest)/(velocity@velocity)
        h=Hole(1,spin=spin,velocity=velocity);sample=lambda x:self.b.seed(h,x)
        out=sample([lab]);self.assertTrue(np.all(np.isfinite(out['Kij'])))
        np.testing.assert_allclose(out['Kij'].reshape(3,3),out['Kij'].reshape(3,3).T,atol=1e-14)
        seq=[norms(constraints(sample,[lab],step)) for step in (.004,.002,.001)]
        self.assertLess(seq[-1]['H_rms'],1e-7);self.assertLess(seq[-1]['M_rms'],1e-7)
        self.assertGreater(seq[0]['M_rms']/seq[1]['M_rms'],10)

if __name__=='__main__':unittest.main()
