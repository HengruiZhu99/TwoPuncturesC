"""Identical adaptive verifier points/steps, without a native context."""
from pathlib import Path
import sys
import unittest
import hashlib
import tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/p) for p in ('python','examples','validation')]
from hispid import Config,Hole
from run_validation import points
from calibrate_study_verifier import sampling_plan,verify_saved_plan


class VerifierSamplingTests(unittest.TestCase):
    def configuration(self):
        cfg=Config();cfg.n[:]=[4,8,4]
        speed=float(np.sqrt(.99))
        cfg.hole[0]=Hole(.5,(12.5,0,0),velocity=(-speed,0,0))
        cfg.hole[1]=Hole(.5,(-12.5,0,0),velocity=(speed,0,0));cfg.inner_max[:]=[.0075]*2
        xyz,near,bulk=points(cfg,True)
        steps=np.minimum(.002,.001*np.min([np.linalg.norm(xyz-np.array(h.center),axis=1) for h in cfg.hole],axis=0))
        return cfg,dict(horizon_scaled=True,verifier_steps=steps.tolist()),xyz,near,bulk,steps

    def test_regions_keep_exact_pointwise_steps(self):
        cfg,record,xyz,near,bulk,steps=self.configuration()
        x,roles,sequence=sampling_plan(cfg,record,'exterior',[2.,1.,.5])
        np.testing.assert_array_equal(x,xyz[:near+bulk]);self.assertTrue(np.all(roles<2))
        np.testing.assert_array_equal(sequence[1],steps[:near+bulk])
        x,roles,sequence=sampling_plan(cfg,record,'all',[2.,1.,.5])
        np.testing.assert_array_equal(x,xyz);self.assertTrue(np.any(roles==2))
        self.assertLess(np.min(sequence[1]),.002)

    def test_duplicate_nonpositive_and_mismatched_steps_are_rejected(self):
        cfg,record,*_=self.configuration()
        for factors in ([2.,1.,1.],[2.,1.,0.],[1.,2.,3.],[2.,float('nan'),.5]):
            with self.subTest(factors=factors),self.assertRaises(ValueError):sampling_plan(cfg,record,'all',factors)
        record['verifier_steps'].pop()
        with self.assertRaises(ValueError):sampling_plan(cfg,record,'all',[2.,1.,.5])

    def test_hash_bound_saved_points_and_steps_must_match(self):
        cfg,record,xyz,_,_,steps=self.configuration()
        record.update(case='synthetic',resolution=list(cfg.n))
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'synthetic_4_4.npz'
            def write(x,h):
                np.savez(path,points=x,verifier_steps=h)
                record['raw_artifact_sha256']={str(path):hashlib.sha256(path.read_bytes()).hexdigest()}
            write(xyz,steps);verify_saved_plan(cfg,record,directory)
            shifted=xyz.copy();shifted[0,0]+=.001
            write(shifted,steps)
            with self.assertRaises(ValueError):verify_saved_plan(cfg,record,directory)
            changed=steps.copy();changed[0]*=2
            write(xyz,changed)
            with self.assertRaises(ValueError):verify_saved_plan(cfg,record,directory)


if __name__=='__main__':unittest.main()
