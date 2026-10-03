"""A failing post-solve callback must leave bound diagnostics and coefficients.

Uses a synthetic Python context, without loading or solving a native system.
"""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/p) for p in ('python','examples','validation')]
from hispid import Config,Hole
from run_validation import solve_case
import execution


class SyntheticContext:
    def __init__(self,config):self.size=4*int(np.prod(config.n));self.resolved_options={}
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def solve(self,**kwargs):return dict(status=0,converged=False,scaled_linf=[2e-14]*4)
    def unknowns(self):return np.zeros(self.size)
    def setup_statistics(self):return None
    def linear_history(self):raise ValueError('synthetic metadata callback failure')


class SyntheticBackend:
    lib=object()
    def library_sha256(self):return 'a'*64
    def residual_scaling(self):return 'sin3_alpha_beta'
    def parameterization(self):return 'synthetic_only'
    def parameterization_description(self):return 'synthetic_only'
    def parameterization_maps(self):return dict(radial_stretch=.2,angular_stretch=2.)
    def create(self,config,execution,geometry='host'):return SyntheticContext(config)


class StudyRetentionTests(unittest.TestCase):
    def test_diagnostic_and_coefficients_survive_history_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);report=root/'results.json';raw=root/'raw'
            def factory(backend,n,nphi):
                config=Config();config.n[:]=[n,n,nphi];config.tolerance=1e-14
                config.hole[0]=Hole(.5,(6,0,0));config.hole[1]=Hole(.5,(-6,0,0))
                return config
            with patch.multiple(execution,name=lambda lib:'Serial',concurrency=lambda lib:1,
                    device_description=lambda lib:None,statistics=lambda lib:{}):
                with self.assertRaisesRegex(ValueError,'synthetic metadata'):
                    solve_case(SyntheticBackend(),factory,[(4,4)],'retention',execution='kokkos',
                        output_report=report,raw_directory=raw)
            saved=json.loads(report.read_text())['retention']
            self.assertFalse(saved['passed']);self.assertEqual(saved['records'],[])
            attempt=saved['incomplete_attempt'];self.assertEqual(attempt['stage'],'collecting_solver_metadata')
            self.assertFalse(attempt['record']['diagnostics']['converged'])
            artifact=attempt['record']['solve_artifact'];path=Path(artifact['path'])
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),artifact['sha256'])
            with np.load(path) as data:self.assertTrue(np.all(data['unknowns']==0))


if __name__=='__main__':unittest.main()
