"""Reject incomplete, nonfinite and amplified ADM refinement evidence.

These tests use analytic arrays only and never construct a native context.
"""
import copy
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'validation'))
from charge_checks import validate_refinements,radial_fits,angular_changes,qualify_refinement


class ChargeRefinementTests(unittest.TestCase):
    radii=[100.,200.,400.,800.,1600.]

    def evidence(self):
        records=[]
        x=1/np.asarray(self.radii)
        values=np.arange(7)[None,:]+2*x[:,None]+3*x[:,None]**2
        for nt,np_ in ((20,40),(32,40),(32,64)):
            item=dict(ntheta=nt,nphi=np_,native_EPJ=values.tolist(),independent_EPJ=values.tolist(),
                native_vs_independent_linf=0.,centered_rotation_linf=0.,
                native_vs_independent_radial_intercepts_linf=0.,fixed_origin_extrapolation_linf=0.)
            radial_fits(item,self.radii,values);radial_fits(item,self.radii,values,'native_')
            angular_changes(item,records[-1] if records else None);records.append(item)
        return records

    def test_separate_refinement_is_required(self):
        validate_refinements(self.radii,[(20,40),(32,40),(32,64)],1e-5)
        for quads in ([(20,40)]*3,[(20,40),(32,64),(48,96)],[(32,64),(20,64),(20,96)]):
            with self.subTest(quads=quads),self.assertRaises(ValueError):
                validate_refinements(self.radii,quads,1e-5)
        with self.assertRaises(ValueError):validate_refinements(self.radii[:4],[(20,40),(32,40),(32,64)],1e-5)

    def test_missing_nonfinite_and_origin_checks_cannot_qualify(self):
        records=self.evidence()
        self.assertTrue(qualify_refinement(records,self.radii,1e-5)['qualified_charge_checks'])
        for key,value in (('fixed_origin_extrapolation_linf',1e-4),
                          ('native_vs_independent_linf',float('nan')),
                          ('radial_fit_windows_EPJ_angular_change_linf',float('inf'))):
            damaged=copy.deepcopy(records);damaged[-1][key]=value
            with self.subTest(key=key):
                self.assertFalse(qualify_refinement(damaged,self.radii,1e-5)['qualified_charge_checks'])
        del records[-1]['native_radial_fit_change_EPJ_linf']
        self.assertFalse(qualify_refinement(records,self.radii,1e-5)['qualified_charge_checks'])

    def test_small_raw_error_can_fail_after_extrapolation(self):
        records=self.evidence();changed=copy.deepcopy(records[-1])
        values=np.asarray(changed['independent_EPJ'])
        # A bounded perturbation at finite radii is amplified by the fit.
        weights=np.polynomial.polynomial.polyfit(1/np.asarray(self.radii[-4:]),np.eye(4),2)[0]
        values[-4:,0]+=9e-6*np.sign(weights)
        changed['independent_EPJ']=values.tolist();radial_fits(changed,self.radii,values)
        angular_changes(changed,records[-2]);records[-1]=changed
        self.assertLess(changed['independent_EPJ_angular_change_linf'],1e-5)
        self.assertGreater(changed['radial_fit_windows_EPJ_angular_change_linf'],1e-5)
        self.assertFalse(qualify_refinement(records,self.radii,1e-5)['qualified_charge_checks'])

    def test_later_joint_step_and_forged_direction_flags_cannot_qualify(self):
        records=self.evidence()
        for record in records:
            record['refined_polar_only']=False;record['refined_azimuthal_only']=False
        self.assertTrue(qualify_refinement(records,self.radii,1e-5)['qualified_charge_checks'])
        last=copy.deepcopy(records[-1]);last.update(ntheta=48,nphi=96)
        last['independent_EPJ_angular_change_linf']=1e-4
        records.append(last)
        self.assertFalse(qualify_refinement(records,self.radii,1e-5)['qualified_charge_checks'])

    def test_reported_all_radius_intercept_must_match_qualified_window(self):
        self.radii=[50.,100.,200.,400.,800.,1600.]
        records=self.evidence()
        for i,record in enumerate(records):
            for key,prefix in (('native_EPJ','native_'),('independent_EPJ','')):
                values=np.asarray(record[key]);values[0,0]+=1e-3
                record[key]=values.tolist();radial_fits(record,self.radii,values,prefix)
            angular_changes(record,records[i-1] if i else None)
        self.assertLess(records[-1]['radial_fit_change_EPJ_linf'],1e-5)
        self.assertGreater(records[-1]['radial_reported_fit_vs_finest_window_linf'],1e-5)
        self.assertFalse(qualify_refinement(records,self.radii,1e-5)['qualified_charge_checks'])


if __name__=='__main__':unittest.main()
