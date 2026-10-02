"""Independent physical-coordinate controls for changed-map initial guesses."""
import unittest
import numpy as np
from prolong import remap_modal,prolong_modal


def physical_values(shape,maps):
    a,b,p=shape;radial,angular=maps
    sigma=.5*(1-np.cos(np.pi*(np.arange(a)+.5)/a))
    t=radial*sigma/(1-(1-radial)*sigma)
    eta=np.tanh(angular*-np.cos(np.pi*(np.arange(b)+.5)/b))/np.tanh(angular)
    v=np.zeros((p,b,a,4))
    polynomial=.003*(1+t[None,:]*eta[:,None]+t[None,:]**2*eta[:,None]**2)
    for k in (0,1,p//2,p//2+1):v[k]=polynomial[:,:,None]*np.array([1.,.4,-.2,.3])
    return v.ravel()


class ModalRemapTests(unittest.TestCase):
    def test_changed_maps_converge_to_physical_polynomial(self):
        oldmaps=(.2,2.);newmaps=(.35,3.);errors=[]
        for n in (16,32,64):
            old=(n,2*n,8);new=(n+5,2*n+7,8)
            remapped=remap_modal(physical_values(old,oldmaps),old,new,oldmaps,newmaps)
            errors.append(float(np.max(abs(remapped-physical_values(new,newmaps)))))
        print('changed-map physical polynomial errors',errors,flush=True)
        self.assertLess(errors[-1],1e-12)
        self.assertLess(errors[1],errors[0]);self.assertLess(errors[2],errors[0])

    def test_same_maps_match_existing_prolongation(self):
        old=(8,10,8);new=(13,17,12);maps=(.2,2.)
        values=physical_values(old,maps)
        np.testing.assert_array_equal(remap_modal(values,old,new,maps,maps),prolong_modal(values,old,new))

    def test_zero_and_nyquist_fourier_normalization(self):
        old=(8,10,8);new=(13,17,12);v=np.zeros((8,10,8,4));v[0,:,:,0]=.03;v[4,:,:,0]=.02;v[5,:,:,1]=.04
        mapped=remap_modal(v.ravel(),old,new,(.2,2.),(.2,3.)).reshape(12,17,13,4)
        phi=2*np.pi*(np.arange(23)+.37)/23
        expected=.03/np.sqrt(8)+.02/np.sqrt(8)*np.cos(4*phi)
        actual=mapped[0,0,0,0]/np.sqrt(12)+mapped[4,0,0,0]*np.sqrt(2/12)*np.cos(4*phi)
        np.testing.assert_allclose(actual,expected,atol=1e-15,rtol=1e-14)
        np.testing.assert_allclose(mapped[7,:,:,1]*np.sqrt(2/12),.04*np.sqrt(2/8),atol=1e-15)
        np.testing.assert_array_equal(mapped[6],0)

    def test_invalid_or_discarded_maps_modes_rejected(self):
        for maps in ((0,2),(.2,0),(.2,float('nan'))):
            with self.assertRaises(ValueError):remap_modal(np.zeros(8*8*8*4),(8,8,8),(8,8,8),maps,maps)
        with self.assertRaises(ValueError):remap_modal(np.zeros(8*8*8*4),(8,8,8),(8,8,4),(.2,2),(.2,3))
        valid=np.zeros(8*8*8*4)
        for shape in ((8,8),(8.,8,8),(0,8,8),(8,8,7),(True,8,8)):
            with self.assertRaises(ValueError):remap_modal(valid,shape,(8,8,8),(.2,2),(.2,2))
            with self.assertRaises(ValueError):remap_modal(valid,(8,8,8),shape,(.2,2),(.2,3))
        for values in (valid[:-1],np.full(valid.shape,np.nan),np.full(valid.shape,np.inf)):
            with self.assertRaises(ValueError):remap_modal(values,(8,8,8),(8,8,8),(.2,2),(.2,2))

if __name__=='__main__':unittest.main()
