"""Prerequisite mutation/compilation guards, without loading native images."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/name) for name in ('python','examples','validation')]
from run_extreme_kokkos import frozen_inputs,prerequisites,verify_hashes


class ExtremePrerequisiteTests(unittest.TestCase):
    def test_frozen_decode_rejects_later_input_mutations(self):
        with tempfile.TemporaryDirectory() as directory:
            paths={key:Path(directory)/key for key in ('plan','performance','seed','receipt','report')}
            for key,path in paths.items():path.write_text('standalone source' if key=='report' else json.dumps({'key':key}))
            payloads,hashes=frozen_inputs(paths)
            self.assertEqual(payloads['plan'],{'key':'plan'})
            verify_hashes(paths,hashes)
            for key,path in paths.items():
                before=path.read_bytes();path.write_bytes(before+b' ')
                with self.subTest(key=key),self.assertRaises(ValueError):verify_hashes(paths,hashes)
                path.write_bytes(before)

    def test_only_complete_matrix_and_exact_compilation_receipt_pass(self):
        hashes={'performance':hashlib.sha256(b'performance').hexdigest(),'report':hashlib.sha256(b'report').hexdigest()}
        performance=dict(declared_performance_completed=True,expected_workers=288,completed_workers=288,
                         records={str(i):{} for i in range(287)},failures={'287':{}})
        receipt=dict(compilation_confirmed=True,compiler='mcp__codex_app__compile_latex_document',
                     performance_sha256=hashes['performance'],report_sha256=hashes['report'])
        prerequisites(performance,receipt,hashes)  # Retained failures do not become successful gates.
        for key,value in (('declared_performance_completed',False),('completed_workers',287),('expected_workers',287)):
            with self.subTest(key=key),self.assertRaises(ValueError):prerequisites(performance|{key:value},receipt,hashes)
        for key,value in (('compilation_confirmed',False),('compiler','unverified'),('performance_sha256','0'*64),('report_sha256','0'*64)):
            with self.subTest(key=key),self.assertRaises(ValueError):prerequisites(performance,receipt|{key:value},hashes)


if __name__=='__main__':unittest.main()
