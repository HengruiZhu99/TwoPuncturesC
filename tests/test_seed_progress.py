"""File/metadata retention controls; no native image or physical calculation."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from check_target_seeds import write_progress


class SeedProgressTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.root=Path(self.directory.name)
        self.source=self.root/'source.py';self.source.write_text('bound implementation\n')
        self.bound={str(self.source):hashlib.sha256(self.source.read_bytes()).hexdigest()}
        self.path=self.root/'attempt.json'
        self.path.write_text('{"passed":false,"completed":false}\n')

    def tearDown(self):self.directory.cleanup()

    def test_partial_update_and_failure_do_not_discard_measurements(self):
        output=dict(passed=False,completed=False,cases=[dict(case='example',passed=False,
            completed=False,charge_quadratures=[dict(charges=[[1,2,3]],completed=False)])])
        metadata=dict(unknown_parameterization_id='bound_basis')
        write_progress(self.path,self.bound,metadata,output)
        output['failure']=dict(type='InterruptedError',message='synthetic interruption')
        write_progress(self.path,self.bound,metadata,output)
        actual=json.loads(self.path.read_bytes())
        self.assertFalse(actual['passed']);self.assertFalse(actual['completed'])
        self.assertEqual(actual['cases'],output['cases'])
        self.assertEqual(actual['failure'],output['failure'])
        self.assertEqual(actual['bound_artifacts_sha256'],self.bound)
        self.assertEqual(actual['unknown_parameterization_id'],'bound_basis')
        self.assertFalse(self.path.with_name('attempt.json.tmp').exists())

    def test_changed_input_preserves_last_bound_snapshot(self):
        previous=self.path.read_bytes();self.source.write_text('changed implementation\n')
        with self.assertRaisesRegex(ValueError,'bound exact-seed'):
            write_progress(self.path,self.bound,{},dict(passed=True,completed=True))
        self.assertEqual(self.path.read_bytes(),previous)

    def test_existing_temporary_does_not_overwrite_evidence(self):
        temporary=self.path.with_name('attempt.json.tmp');temporary.write_text('retained partial write')
        previous=self.path.read_bytes()
        with self.assertRaises(FileExistsError):
            write_progress(self.path,self.bound,{},dict(passed=True,completed=True))
        self.assertEqual(self.path.read_bytes(),previous)
        self.assertEqual(temporary.read_text(),'retained partial write')

    def invoke_missing_library(self,path):
        script=Path(write_progress.__code__.co_filename).resolve()
        return subprocess.run([sys.executable,str(script),'--library',str(self.root/'missing.so'),
            '--output',str(path)],capture_output=True,text=True)

    def test_setup_failure_has_explicit_unqualified_record(self):
        path=self.root/'setup.json';result=self.invoke_missing_library(path)
        self.assertNotEqual(result.returncode,0)
        actual=json.loads(path.read_bytes())
        self.assertFalse(actual['passed']);self.assertFalse(actual['completed'])
        self.assertEqual(actual['stage'],'setup_failed')
        self.assertEqual(actual['failure']['type'],'FileNotFoundError')

    def test_unavailable_failure_write_preserves_original_setup_error(self):
        path=self.root/'setup.json';temporary=path.with_name(path.name+'.tmp')
        temporary.write_text('retained temporary evidence')
        result=self.invoke_missing_library(path)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('FileNotFoundError',result.stderr.splitlines()[-1])
        self.assertIn('prior snapshot preserved',result.stderr)
        actual=json.loads(path.read_bytes())
        self.assertFalse(actual['passed']);self.assertFalse(actual['completed'])
        self.assertEqual(actual['stage'],'loading_bound_library')
        self.assertEqual(temporary.read_text(),'retained temporary evidence')


if __name__=='__main__':unittest.main()
