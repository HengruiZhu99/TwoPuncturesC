"""Compiled-source and dependency witnesses for the separate setup manifest."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/name) for name in ('python','examples','validation')]
from kokkos_build_manifest import bind_build,native_source,setup_build_source
from benchmark_bowen_york import digest
from perlmutter_allocation_guard import expected_workers


class SetupBuildManifestTests(unittest.TestCase):
    def test_shared_guard_keeps_both_protocol_inventories(self):
        legacy=dict(binding={'grids':[1,2,3],'repeats':3,'systems':['hispid','by'],
                             'methods':['gmres','bicgstab']},manifest={'variants':list(range(8))})
        self.assertEqual(expected_workers(legacy),288)
        setup=dict(expected_worker_ids=['host','execution'],records={'host':{}},failures={})
        self.assertEqual(expected_workers(setup),2)
        for changed in (setup|{'expected_worker_ids':['host','host']},
                        setup|{'expected_worker_ids':[]}, setup|{'records':{'unexpected':{}}},
                        setup|{'failures':{'host':{}}}):
            with self.assertRaises(ValueError):expected_workers(changed)

    def fixture(self,directory):
        root=Path(directory);source=root/'source';build=root/'build';source.mkdir();build.mkdir()
        for name in ('src','include','tests'):(source/name).mkdir()
        (source/'CMakeLists.txt').write_text('project(control)\n')
        (source/'src/native.cpp').write_text('int control(){return 1;}\n')
        cache='# This is the CMakeCache file\n\n// A cache entry\n'+'\n\n// A cache entry\n'.join(('CMAKE_BUILD_TYPE:STRING=Release','PUNCTURES_KOKKOS:BOOL=ON',
                         'PUNCTURES_BENCHMARK:BOOL=ON','HISPID_ROW_POWER:STRING=3',
                         'CMAKE_C_FLAGS:STRING=','CMAKE_CXX_FLAGS:STRING=',
                         'CMAKE_C_FLAGS_RELEASE:STRING=-O3 -DNDEBUG',
                         'CMAKE_CXX_FLAGS_RELEASE:STRING=-O3 -DNDEBUG',
                         'CMAKE_HOME_DIRECTORY:INTERNAL='+str(source)))+'\n'
        (build/'CMakeCache.txt').write_text(cache)
        for name in ('libHiSpID.so','libTwoPunctures.so','libkokkoscore.so'):
            (build/name).write_bytes(name.encode())
        listing='libkokkoscore.so => '+str(build/'libkokkoscore.so')+' (0x0000)\n'
        return source,build,cache,listing

    def test_bind_actual_dependencies_and_native_source(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'validation') as directory:
            source,build,cache,listing=self.fixture(directory)
            with patch('kokkos_build_manifest.subprocess.check_output',return_value=listing):
                record,runtime,dependencies=bind_build(build)
            self.assertEqual(runtime,[str(build/'libkokkoscore.so')])
            self.assertEqual(record['images'][str(build/'libHiSpID.so')],digest(build/'libHiSpID.so'))
            self.assertEqual(record['cache_sha256'],digest(build/'CMakeCache.txt'))
            expected=native_source(source)
            self.assertEqual(setup_build_source(record,expected)['native_source_sha256'],expected)
            (source/'src/native.cpp').write_text('int control(){return 2;}\n')
            with self.assertRaisesRegex(ValueError,'compiled native source differs'):
                setup_build_source(record,expected)

    def test_reject_wrong_flags_missing_or_unresolved_dependencies(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'validation') as directory:
            source,build,cache,listing=self.fixture(directory)
            for replacement in ('HISPID_ROW_POWER:STRING=6','PUNCTURES_KOKKOS:BOOL=OFF',
                                'CMAKE_BUILD_TYPE:STRING=Debug','PUNCTURES_BENCHMARK:BOOL=OFF'):
                key=replacement.split(':')[0]
                altered='\n'.join(replacement if line.startswith(key+':') else line for line in cache.splitlines())
                with self.subTest(replacement=replacement),self.assertRaises(ValueError):
                    setup_build_source({'cache':altered},native_source(source))
            for flag in ('-ffast-math','-Ofast','--use_fast_math'):
                with self.subTest(flag=flag),self.assertRaisesRegex(ValueError,'fast math'):
                    setup_build_source({'cache':cache.replace('-O3 -DNDEBUG',flag)},native_source(source))
            with patch('kokkos_build_manifest.subprocess.check_output',return_value='libkokkoscore.so => not found\n'):
                with self.assertRaisesRegex(RuntimeError,'unresolved dependency'):
                    bind_build(build)
            (build/'libHiSpID.so').unlink()
            with patch('kokkos_build_manifest.subprocess.check_output',return_value=listing):
                with self.assertRaises(FileNotFoundError):bind_build(build)


if __name__=='__main__':unittest.main()
