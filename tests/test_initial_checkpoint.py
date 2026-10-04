"""Reject reinterpretation of portable guesses before creating native contexts."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/name) for name in ('python','examples','validation')]
from hispid import Config
from checkpoint_export import write_checkpoint
from run_extreme_kokkos import binary_config
from run_validation import solve_case


class ConfigOnlyBackend:
    context_calls=0
    def config(self):return Config()
    def library_sha256(self):return 'a'*64
    def parameterization(self):return 'modal_P_C2prolate_mapped_v2'
    def parameterization_description(self):return 'guard-only dummy'
    def parameterization_maps(self):return dict(radial_stretch=.2,angular_stretch=2.)
    def residual_scaling(self):return 'sin3_alpha_beta'
    def create(self,*args,**kwargs):
        self.context_calls+=1
        raise RuntimeError('valid guess reached context creation')


class PortableInitialGuessGuards(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory();self.root=Path(self.directory.name)
        self.backend=ConfigOnlyBackend();self.path=self.root/'guess.checkpoint'
        self.controls=dict(max_newton=24,max_krylov=4800,memory_limit_mib=32768,
                           outer_tolerance=1e-14,restart=200)
        self.cfg=binary_config(self.backend,[4,4,4],'headon_gamma10_kokkos',25.,self.controls)
        self.values=np.linspace(-.1,.1,256)
        self.write()

    def tearDown(self):self.directory.cleanup()

    def write(self,cfg=None,producer='a'*64,token='modal_P_C2prolate_mapped_v2'):
        write_checkpoint(self.path,self.cfg if cfg is None else cfg,self.values,producer,'diagnostic',token)

    def run_case(self,**options):
        return solve_case(self.backend,lambda backend,n,nphi:self.cfg,[(4,4)],'fresh',
                          output_report=self.root/'report.json',raw_directory=self.root/'raw',
                          initial_checkpoint=self.path,**options)

    def assert_rejected(self):
        with self.assertRaises(ValueError):self.run_case()
        self.assertEqual(self.backend.context_calls,0)

    def test_wrong_producer_or_basis_never_creates_context(self):
        self.write(producer='b'*64);self.assert_rejected()
        self.write(token='W_plus_Aminus1_V');self.assert_rejected()

    def test_changed_boost_or_grid_never_creates_context(self):
        different=Config.from_buffer_copy(bytes(self.cfg));different.hole[0].velocity[0]=-.5
        self.write(different);self.assert_rejected()
        different=Config.from_buffer_copy(bytes(self.cfg));different.n[:]=[4,8,4]
        self.values=np.zeros(512);self.write(different);self.assert_rejected()

    def test_existing_result_and_other_warm_start_are_rejected(self):
        (self.root/'report.json').write_text(json.dumps({'fresh':{}}));self.assert_rejected()
        (self.root/'report.json').unlink()
        with self.assertRaises(ValueError):self.run_case(initial_guess='other-guess.json')
        self.assertEqual(self.backend.context_calls,0)

    def test_solver_control_changes_are_allowed_only_for_a_fresh_context(self):
        source=Config.from_buffer_copy(bytes(self.cfg));source.krylov_restart=80;source.max_krylov=2400
        self.write(source)
        with self.assertRaisesRegex(RuntimeError,'valid guess reached context creation'):self.run_case()
        self.assertEqual(self.backend.context_calls,1)


if __name__=='__main__':unittest.main()
