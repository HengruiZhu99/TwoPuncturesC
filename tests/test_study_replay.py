"""Separate-study provenance and covariance preflight, without native contexts."""
import copy
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/p) for p in ('python','examples','validation')]
from hispid import Config
from configs import as_dict
from checkpoints import restore_payload,PARAMETERIZATION
from check_covariance import validate_source_sequence


class StudyReplayTests(unittest.TestCase):
    def test_relocated_raw_bytes_are_bound_and_tampering_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            raw=Path(directory);cfg=Config();cfg.n[:]=[4,8,4]
            record=dict(case='separate',resolution=list(cfg.n),config=as_dict(cfg),unknown_parameterization=PARAMETERIZATION)
            path=raw/'separate_4_4.npz';unknowns=np.arange(4*4*8*4,dtype=float)
            np.savez_compressed(path,unknowns=unknowns)
            record['raw_artifact_sha256']={'/old/location/'+path.name:hashlib.sha256(path.read_bytes()).hexdigest()}
            restored,values=restore_payload(record,Config(),raw)
            self.assertEqual(list(restored.n),list(cfg.n));np.testing.assert_array_equal(values,unknowns)
            np.savez_compressed(path,unknowns=unknowns+1)
            with self.assertRaisesRegex(ValueError,'retained evidence'):restore_payload(record,Config(),raw)

    def test_full_grid_free_data_and_sampling_must_match(self):
        levels=[(40,16),(80,20),(128,28)]
        records=[dict(resolution=[n,2*n,p],config=dict(n=[n,2*n,p],hole=[{'spin':[0,0,.99]}]),horizon_scaled=True)
                 for n,p in levels]
        self.assertEqual(validate_source_sequence(records,levels),records)
        mutations=[]
        polar=copy.deepcopy(records);polar[1]['resolution'][1]=64;polar[1]['config']['n'][1]=64;mutations.append(polar)
        spin=copy.deepcopy(records);spin[1]['config']['hole'][0]['spin'][2]=.95;mutations.append(spin)
        points=copy.deepcopy(records);points[1]['horizon_scaled']=False;mutations.append(points)
        for rows in mutations:
            with self.subTest(rows=rows),self.assertRaises(ValueError):validate_source_sequence(rows,levels)


if __name__=='__main__':unittest.main()
