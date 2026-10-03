"""Reject invalid supplemental setup data without executing native kernels."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/name) for name in ('python', 'examples', 'validation')]
from benchmark_geometry_setup import SETUP_TESTS, compare_arrays, decode_config, finalize, frozen_arrays, junit_outcomes, protocol_checks, state_shapes, verify_sources, verify_controls
from configs import as_dict
from hispid import Config, Hole


class GeometrySetupProtocolTests(unittest.TestCase):
    def fixtures(self, geometry='execution', space='Serial'):
        variant = dict(id='serial', threads=1, space=space, execution='kokkos',
                       hispid_library='/build/libHiSpID.so', by_library='/build/libTwoPunctures.so', runtime_images=['/build/libkokkoscore.so'])
        images = {'/build/libHiSpID.so': 'a'*64, '/build/libTwoPunctures.so': 'b'*64,
                  '/build/libkokkoscore.so': 'c'*64}
        record = dict(config={'n': [4, 8, 4]}, geometry=geometry, execution='kokkos',
                      input_sha256='d'*64, input_arrays_sha256='e'*64,
                      compiled_execution=space, cpu_threads=1,
                      library_sha256='a'*64, loaded_image_verified=True,
                      dependency_images={'/build/libTwoPunctures.so': 'b'*64},
                      runtime_images={'/build/libkokkoscore.so': 'c'*64},
                      all_arrays_finite=True, initial_unknowns_zero=True,
                      diagnostics={'newton_iterations': 0, 'krylov_iterations': 0},
                      setup_statistics={'geometry_execution': int(geometry == 'execution'),
                                        'scalar_digits': 53 if geometry == 'execution' else 64,
                                        'spectral_seconds': .1, 'geometry_seconds': .2, 'coefficient_seconds': .01},
                      execution_statistics={'memory_tracking_available': 1},
                      initialization_seconds=.01, creation_seconds=.4, first_sample_seconds=.02,
                      ready_to_sample_seconds=.42, driver_memory={'samples': []},
                      device={'visible_count': 1, 'visible_ordinal': 0, 'uuid': 'GPU-example'})
        return variant, images, record

    def checks(self, record, variant, images):
        return protocol_checks(record, variant, images, {'n': [4, 8, 4]}, 'd'*64, 'e'*64)

    def test_frozen_config_and_arrays(self):
        config = Config()
        config.hole[0] = Hole(.6, (3, 0, 0), (.072, .054, .108), (.03, .06, .01))
        config.hole[1] = Hole(.4, (-3, 0, 0), (-.032, .048, .016), (-.02, -.07, .025))
        config.n[:] = [4, 8, 4]
        values = as_dict(config)
        self.assertEqual(as_dict(decode_config(values)), values)
        with self.assertRaises(ValueError):
            decode_config(values | {'extra': 1})
        first, second = frozen_arrays([4, 8, 4]), frozen_arrays([4, 8, 4])
        for key in first:
            np.testing.assert_array_equal(first[key], second[key])
        self.assertEqual(first['unknowns'].shape, (512,))
        self.assertFalse(np.array_equal(first['unknowns'], first['direction']))

    def test_reject_nonfinite_shapes_and_keys(self):
        with tempfile.TemporaryDirectory() as directory:
            host, execution = Path(directory)/'host.npz', Path(directory)/'execution.npz'
            arrays = {key: np.zeros(shape, dtype=float) for key, shape in state_shapes([4,8,4]).items()}
            np.savez(host, **arrays); np.savez(execution, **arrays)
            self.assertTrue(compare_arrays(host, execution, [4,8,4])['passed'])
            for mutation in (arrays | {'zero_residual': np.full(512, np.nan)},
                             arrays | {'zero_residual': np.zeros(4)},
                             arrays | {'zero_residual': np.zeros(512, dtype=np.float32)},
                             arrays | {'zero_residual': np.full(512, 1e-7)},
                             {key: value for key, value in arrays.items() if key != 'dgamma'}):
                np.savez(execution, **mutation)
                self.assertFalse(compare_arrays(host, execution, [4,8,4])['passed'])
            np.savez(host, residual=np.zeros(4)); np.savez(execution, residual=np.zeros(4))
            self.assertFalse(compare_arrays(host, execution, [4,8,4])['passed'])

    def test_protocol_mutations_remain_failed(self):
        variant, images, record = self.fixtures()
        self.assertTrue(self.checks(record, variant, images)['passed'])
        mutations = [dict(config={'n': [8, 8, 4]}), dict(input_arrays_sha256='f'*64),
                     dict(cpu_threads=2), dict(compiled_execution='OpenMP'),
                     dict(library_sha256='f'*64), dict(runtime_images={}), dict(dependency_images={}),
                     dict(all_arrays_finite=False), dict(initial_unknowns_zero=False),
                     dict(diagnostics={'newton_iterations': 1, 'krylov_iterations': 0}),
                     dict(creation_seconds=float('nan')), dict(creation_seconds=0), dict(first_sample_seconds=-1),
                     dict(execution_statistics={'memory_tracking_available': 0}),
                     dict(setup_statistics=record['setup_statistics'] | {'geometry_execution': 0}),
                     dict(setup_statistics=record['setup_statistics'] | {'scalar_digits': 64})]
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                self.assertFalse(self.checks(record | mutation, variant, images)['passed'])

    def test_gpu_identity(self):
        variant, images, record = self.fixtures(space='Cuda')
        self.assertTrue(self.checks(record, variant, images)['passed'])
        for mutation in (dict(device=record['device'] | {'visible_count': 2}),
                         dict(device=record['device'] | {'uuid': None}),
                         dict(driver_memory={'samples': [{'uuid': 'GPU-other'}]})):
            self.assertFalse(self.checks(record | mutation, variant, images)['passed'])

    def test_pair_needs_both_suites_and_both_worker_checks(self):
        variant, images, execution = self.fixtures()
        _, _, host = self.fixtures(geometry='host')
        with tempfile.TemporaryDirectory(dir=ROOT/'validation') as directory:
            raw = Path(directory)
            state = raw/'state.npz'
            arrays = {key: np.zeros(shape) for key, shape in state_shapes([4,8,4]).items()}
            np.savez(state, **arrays)
            import hashlib
            sha = hashlib.sha256(state.read_bytes()).hexdigest()
            for record in (host, execution):
                record.update(variant='serial', grid=[4, 8, 4], repeat=0,
                              state=str(state.relative_to(ROOT)), state_sha256=sha)
            result = dict(records={'4_8_4_serial_0_host': host, '4_8_4_serial_0_execution': execution}, failures={},
                          setup_controls={'serial': dict(passed=True, executed=True, exit_code=0, timeout=False,
                              tests={key: True for key in SETUP_TESTS}, images=images,
                              executables_sha256={'/test/setup': 'f'*64, '/test/seed': 'f'*64},
                              log='control.log', log_sha256='0'*64, junit='control.xml', junit_sha256='0'*64)}, manifest={'variants': [variant]},
                          images={'serial': images}, binding={'grids': [[4, 8, 4]], 'repeats': 1, 'sources_sha256': {}},
                          inputs={'4_8_4': dict(config={'n': [4, 8, 4]}, config_sha256='d'*64, arrays_sha256='e'*64)},
                          expected_worker_ids=['4_8_4_serial_0_host', '4_8_4_serial_0_execution'])
            finalize(result, raw, verify=False)
            self.assertFalse(result['full_artifact_verification'])
            self.assertFalse(result['comparisons']['4_8_4_serial_0']['qualified'])
            self.assertFalse(result['all_setup_checks_passed'])
            with patch('benchmark_geometry_setup.verify_controls'):
                finalize(result, raw)
            self.assertTrue(result['comparisons']['4_8_4_serial_0']['qualified'])
            self.assertTrue(result['all_setup_checks_passed'])
            self.assertFalse(result['declared_setup_completed'])
            for change in ('native_control', 'skipped_suite', 'missing_suite', 'worker', 'array'):
                failed = copy.deepcopy(result)
                if change == 'native_control':
                    failed['setup_controls']['serial']['passed'] = False
                elif change == 'skipped_suite':
                    failed['setup_controls']['serial']['tests']['hispid_kokkos_setup'] = False
                elif change == 'missing_suite':
                    del failed['setup_controls']['serial']['tests']['hispid_kokkos_setup']
                elif change == 'worker':
                    failed['records']['4_8_4_serial_0_execution']['all_arrays_finite'] = False
                else:
                    other = raw/'bad.npz'
                    np.savez(other, **(arrays | {'zero_residual': np.full(512, 1e-7)}))
                    failed['records']['4_8_4_serial_0_execution'].update(state=str(other.relative_to(ROOT)),
                            state_sha256=hashlib.sha256(other.read_bytes()).hexdigest())
                with patch('benchmark_geometry_setup.verify_controls'):
                    finalize(failed, raw)
                comparison = failed['comparisons']['4_8_4_serial_0']
                self.assertFalse(comparison['qualified'])
                self.assertNotIn('creation_speedup', comparison)
                self.assertFalse(failed['all_setup_checks_passed'])
            # A fast last-worker checkpoint must not look qualified if the
            # following full source/control/artifact verification rejects it.
            finalize(result, raw, verify=False)
            with patch('benchmark_geometry_setup.verify_controls', side_effect=ValueError('changed control')):
                with self.assertRaisesRegex(ValueError, 'changed control'):
                    finalize(result, raw)
            self.assertFalse(result['all_setup_checks_passed'])
            self.assertFalse(result['declared_setup_completed'])

    def test_source_and_control_witness_mutations(self):
        import hashlib
        with tempfile.TemporaryDirectory(dir=ROOT/'validation') as directory:
            path = Path(directory)/'witness'; path.write_text('original')
            sha = hashlib.sha256(path.read_bytes()).hexdigest()
            result = dict(binding={'sources_sha256': {str(path.relative_to(ROOT)): sha}},
                          setup_controls={'serial': {'executables_sha256': {str(path): sha}, 'images': {str(path): sha}}})
            verify_sources(result); verify_controls(result)
            path.write_text('changed')
            with self.assertRaisesRegex(ValueError, 'source changed'):
                verify_sources(result)
            with self.assertRaisesRegex(ValueError, 'executable changed'):
                verify_controls(result)
            result['setup_controls']['serial']['executables_sha256'] = {}
            with self.assertRaisesRegex(ValueError, 'image changed'):
                verify_controls(result)

    def test_relocated_artifacts_and_fresh_junit_outcomes(self):
        import hashlib
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, log, junit, image = [root/name for name in ('source.py', 'control.log', 'control.xml', 'native')]
            source.write_text('source'); log.write_text('control'); image.write_text('image')
            junit.write_text('<testsuite><testcase name="hispid_kokkos_setup"/><testcase name="hispid_execution_seed_export"/></testsuite>')
            sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
            row = dict(log='control.log', log_sha256=sha(log), junit='control.xml', junit_sha256=sha(junit),
                       tests={key: True for key in SETUP_TESTS}, executables_sha256={str(image): sha(image)},
                       images={str(image): sha(image)})
            result = dict(binding={'sources_sha256': {'source.py': sha(source)}}, setup_controls={'serial': row})
            verify_sources(result, artifact_root=root); verify_controls(result, artifact_root=root)
            junit.write_text('<testsuite><testcase name="hispid_kokkos_setup"><skipped/></testcase><testcase name="hispid_execution_seed_export"/></testsuite>')
            row['junit_sha256'] = sha(junit)
            with self.assertRaisesRegex(ValueError, 'outcomes differ'):
                verify_controls(result, artifact_root=root)
            row['tests']['hispid_kokkos_setup'] = False
            verify_controls(result, artifact_root=root)
            junit.write_text('<testsuite><testcase name="hispid_kokkos_setup"><failure/></testcase><testcase name="hispid_kokkos_setup"/><testcase name="hispid_execution_seed_export"/></testsuite>')
            row['junit_sha256'] = sha(junit)
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                verify_controls(result, artifact_root=root)
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                junit_outcomes(junit)


if __name__ == '__main__':
    unittest.main()
