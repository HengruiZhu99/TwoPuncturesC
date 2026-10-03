"""Report qualification controls using synthetic files, without native runs."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/name) for name in ('python', 'examples', 'validation')]
from benchmark_geometry_setup import BOUND, GRIDS, SETUP_TESTS, frozen_arrays, state_shapes
from configs import as_dict
from hispid import Config
from kokkos_latex_report import number, qualified_solve_group, tex
from kokkos_setup_report import SPACES, inventory, load_setup_results, qualified_group, render_setup, verify_manifest


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


class SetupReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        root = cls.root = Path(cls.temporary.name).resolve()
        image_paths = [root/name for name in ('libHiSpID.so', 'libTwoPunctures.so', 'libkokkoscore.so')]
        for image in image_paths:
            image.write_text('synthetic image: '+image.name)
        images = {str(p): sha(p) for p in image_paths}
        variants = [dict(id=key, space=space, threads=threads, execution='kokkos',
                         hispid_library=str(image_paths[0]), by_library=str(image_paths[1]),
                         runtime_images=[str(image_paths[2])]) for key, (space, threads) in SPACES.items()]
        log, junit = root/'control.log', root/'control.xml'
        log.write_text('synthetic native-control metadata')
        junit.write_text('<testsuite>'+''.join(f'<testcase name="{name}"/>' for name in sorted(SETUP_TESTS))+'</testsuite>')
        control = dict(executed=True, passed=True, tests={key: True for key in SETUP_TESTS}, exit_code=0, timeout=False,
                       images=images, executables_sha256={str(p): sha(p) for p in image_paths[:2]},
                       log=log.name, log_sha256=sha(log), junit=junit.name, junit_sha256=sha(junit))
        result = dict(manifest={'variants': variants}, manifest_sha256='a'*64,
                      binding=dict(grids=GRIDS, repeats=3, timeout=3600, bound=BOUND, sources_sha256={},
                                   geometry_choices=['host', 'execution'], nonlinear_solve=False),
                      expected_worker_ids=[f'{"_".join(map(str,g))}_{v["id"]}_{repeat}_{geometry}' for g in GRIDS
                                           for v in variants for repeat in range(3) for geometry in ('host', 'execution')],
                      records={}, failures={}, images={v['id']: images for v in variants}, setup_controls={'serial': control},
                      inputs={}, input_sha256={}, comparisons={})
        for grid in GRIDS:
            key = '_'.join(map(str, grid))
            config = Config(); config.n[:] = grid
            values = as_dict(config)
            config_path, array_path = root/(key+'.json'), root/(key+'.npz')
            save(config_path, values); np.savez(array_path, **frozen_arrays(grid))
            result['inputs'][key] = dict(config=values, config_path=config_path.name, config_sha256=sha(config_path),
                                        arrays_path=array_path.name, arrays_sha256=sha(array_path))
            result['input_sha256'].update({p.name: sha(p) for p in (config_path, array_path)})
        grid = GRIDS[0]
        cls.arrays = {key: np.zeros(shape) for key, shape in state_shapes(grid).items()}
        state = root/'state.npz'; np.savez(state, **cls.arrays)
        worker_log = root/'worker.log'; worker_log.write_text('synthetic setup worker')
        inputs = result['inputs']['_'.join(map(str,grid))]
        for repeat in range(3):
            for geometry in ('host', 'execution'):
                label = f'{"_".join(map(str,grid))}_serial_{repeat}_{geometry}'
                create = 4. if geometry == 'host' else 2.
                worker = dict(config=inputs['config'], geometry=geometry, execution='kokkos',
                              compiled_execution='Serial', cpu_threads=1, loaded_image_verified=True,
                              library_sha256=images[str(image_paths[0])],
                              dependency_images={str(image_paths[1]): images[str(image_paths[1])]},
                              runtime_images={str(image_paths[2]): images[str(image_paths[2])]},
                              input_sha256=inputs['config_sha256'], input_arrays_sha256=inputs['arrays_sha256'],
                              all_arrays_finite=True, initial_unknowns_zero=True,
                              diagnostics={'newton_iterations': 0, 'krylov_iterations': 0},
                              setup_statistics=dict(geometry_execution=int(geometry=='execution'), scalar_digits=53 if geometry=='execution' else 64,
                                                    spectral_seconds=.1, geometry_seconds=create-.1, coefficient_seconds=.05),
                              execution_statistics=dict(memory_tracking_available=1, kokkos_host_peak_bytes=1024, kokkos_device_peak_bytes=0),
                              initialization_seconds=.1, creation_seconds=create, first_sample_seconds=.1,
                              ready_to_sample_seconds=create+.1, instrumentation_gap_seconds=.01,
                              wall_creation_through_sample_seconds=create+.11, max_rss_bytes=8192 if geometry=='host' else 4096)
                worker_path = root/(label+'.json'); save(worker_path, worker)
                result['records'][label] = dict(worker, variant='serial', repeat=repeat, grid=grid,
                                               state=state.name, state_sha256=sha(state), worker=worker_path.name, worker_sha256=sha(worker_path),
                                               log=worker_log.name, log_sha256=sha(worker_log), driver_memory={'samples': [], 'sampled_process_peak_bytes': None})
        cls.template = result

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def load(self, result):
        path = self.root/'results.json'; save(path, result)
        before = path.read_bytes()
        # Binary build validation has its own mutation tests below. These
        # synthetic image/control files make no native or GPU validation claim.
        with patch('kokkos_setup_report.verify_manifest'):
            loaded, digest = load_setup_results(path, self.root)
        self.assertEqual(path.read_bytes(), before)
        return loaded, digest

    def test_fresh_checks_ignore_cached_flags_and_require_three_pairs(self):
        result = copy.deepcopy(self.template)
        result.update(declared_setup_completed=True, all_setup_checks_passed=True)
        result['comparisons'] = {'forged': {'qualified': True, 'passed': True}}
        loaded, digest = self.load(result)
        self.assertFalse(loaded['declared_setup_completed'])
        self.assertFalse(loaded['all_setup_checks_passed'])
        self.assertNotIn('forged', loaded['comparisons'])
        ratio = qualified_group(loaded, GRIDS[0], 'serial')
        self.assertEqual(ratio['creation_speedup'], 2.)
        self.assertEqual(ratio['rss_saved_percent'], 50.)
        report = render_setup(loaded, digest, tex, number)
        self.assertIn('6/126', report)
        self.assertIn('Setup campaign pending', report)
        self.assertIn('No nonlinear', report.replace('It performs no nonlinear', 'No nonlinear'))
        del result['records']['40_80_16_serial_2_execution']
        loaded, _ = self.load(result)
        self.assertIsNone(qualified_group(loaded, GRIDS[0], 'serial'))

    def test_coverage_identity_and_required_evidence(self):
        for mutation in ('inventory', 'identity', 'overlap', 'missing_log', 'variant', 'control'):
            result = copy.deepcopy(self.template)
            first = next(iter(result['records']))
            if mutation == 'inventory': result['expected_worker_ids'].pop()
            elif mutation == 'identity': result['records'][first]['repeat'] = 2
            elif mutation == 'overlap': result['failures'][first] = result['records'][first]
            elif mutation == 'missing_log': del result['records'][first]['log']
            elif mutation == 'variant': result['manifest']['variants'][0]['threads'] = 2
            else: result['setup_controls'] = {}
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                inventory(result)

    def test_record_and_input_witnesses(self):
        result = copy.deepcopy(self.template)
        first = next(iter(result['records']))
        result['records'][first]['creation_seconds'] = 100.
        with self.assertRaisesRegex(ValueError, 'original worker'):
            self.load(result)
        result = copy.deepcopy(self.template)
        result['inputs']['40_80_16']['config']['mass_floor'] = .123
        with self.assertRaisesRegex(ValueError, 'configuration'):
            self.load(result)
        result = copy.deepcopy(self.template)
        result['input_sha256'] = {}
        with self.assertRaisesRegex(ValueError, 'input hashes'):
            self.load(result)

    def test_cached_array_success_cannot_hide_missing_keys(self):
        result = copy.deepcopy(self.template)
        state = self.root/'incomplete.npz'
        np.savez(state, **{k:v for k,v in self.arrays.items() if k != 'dgamma'})
        row = result['records']['40_80_16_serial_0_execution']
        row.update(state=state.name, state_sha256=sha(state))
        result['comparisons']['40_80_16_serial_0'] = dict(passed=True, qualified=True)
        loaded, _ = self.load(result)
        self.assertFalse(loaded['comparisons']['40_80_16_serial_0']['qualified'])
        self.assertIsNone(qualified_group(loaded, GRIDS[0], 'serial'))

    def test_manifest_sources_caches_and_images_are_rehashed(self):
        root = self.root
        source, cache = root/'source.py', root/'CMakeCache.txt'
        source.write_text('synthetic source'); cache.write_text('synthetic cache')
        result = copy.deepcopy(self.template)
        manifest = result['manifest']
        manifest.update(scope='supplemental paired geometry setup', solve_campaign_replaced=False,
                        source_sha256={'source.py': sha(source)},
                        builds={key: dict(cache_path=str(cache), cache_sha256=sha(cache), cache=cache.read_text(), compiled_source={},
                                          images={p:s for p,s in result['images']['serial'].items() if 'libkokkos' not in p},
                                          runtime_images={p:s for p,s in result['images']['serial'].items() if 'libkokkos' in p})
                                for key in ('serial', 'openmp', 'cuda')})
        def bind():
            result['manifest_sha256'] = hashlib.sha256((json.dumps(manifest, indent=2)+'\n').encode()).hexdigest()
        bind()
        with patch('kokkos_setup_report.native_source', return_value={}), patch('kokkos_setup_report.setup_build_source', return_value={}), \
                patch('kokkos_setup_report.manifest_images', return_value=result['images']['serial']):
            verify_manifest(result, root)
            source.write_text('changed')
            with self.assertRaisesRegex(ValueError, 'source changed'): verify_manifest(result, root)
            source.write_text('synthetic source')
            cache.write_text('changed')
            with self.assertRaisesRegex(ValueError, 'cache changed'): verify_manifest(result, root)
            cache.write_text('synthetic cache')
            saved = manifest['builds'].pop('cuda'); bind()
            with self.assertRaisesRegex(ValueError, 'build inventory'): verify_manifest(result, root)
            manifest['builds']['cuda'] = saved; bind()
            saved_images = manifest['builds']['serial']['images']
            manifest['builds']['serial']['images'] = {}; bind()
            with self.assertRaisesRegex(ValueError, 'image bank'): verify_manifest(result, root)
            manifest['builds']['serial']['images'] = saved_images; bind()
            result['images']['serial'] = {}
            with self.assertRaisesRegex(ValueError, 'images changed'): verify_manifest(result, root)
        result['manifest_sha256'] = 'f'*64
        with self.assertRaisesRegex(ValueError, 'manifest differs'): verify_manifest(result, root)

    def test_solve_ratios_require_complete_preserved_repeats(self):
        rows = [{'repeat': r, 'checks': {'passed': True}} for r in range(3)]
        result = dict(binding={'repeats': 3}, comparisons={f'40_80_16_cuda_by_gmres_{r}': {'passed': True} for r in range(3)})
        valid = lambda: qualified_solve_group(result, GRIDS[0], 'by', 'gmres', 'cuda', rows)
        self.assertTrue(valid())
        rows.pop(); self.assertFalse(valid())
        rows.append({'repeat': 2, 'checks': {'passed': True}})
        result['comparisons']['40_80_16_cuda_by_gmres_2']['passed'] = False
        self.assertFalse(valid())
        result['comparisons']['40_80_16_cuda_by_gmres_2']['passed'] = True
        rows[0]['checks']['passed'] = False; self.assertFalse(valid())


if __name__ == '__main__':
    unittest.main()
