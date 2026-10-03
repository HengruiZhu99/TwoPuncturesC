"""Coverage and retained-artifact guards must not accept partial results."""
import copy,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import benchmark_kokkos as b

class ProtocolTests(unittest.TestCase):
    def fixture(self):
        source=json.loads((b.ROOT/'validation/kokkos_preflight40_20261002.json').read_text())
        label='40_80_16_reference_hispid_gmres_0';row=copy.deepcopy(source['records'][label])
        for key in ('state','state_sha256','log','log_sha256'):row.pop(key,None)
        row['runtime_images']={}
        row['resolved_options']=dict(krylov='gmres',linear_rtol=.001,preconditioner='modal',native_verified=False)
        variant=source['manifest']['variants'][0]
        return label,dict(records={label:row},failures={},manifest=dict(variants=[variant]),images=source['images'],configurations={'40:80:16':copy.deepcopy(row['config'])})
    def test_complete(self):
        label,result=self.fixture()
        with patch.object(b,'hi_compare',return_value=dict(passed=True)):b.finalize(result,[label])
        self.assertTrue(result['all_stopping_checks_passed']);self.assertTrue(result['strict_state_comparisons_passed'])
    def test_incomplete(self):
        label,result=self.fixture()
        with patch.object(b,'hi_compare',return_value=dict(passed=True)):b.finalize(result,[label,label[:-1]+'1'])
        self.assertFalse(result['coverage_complete']);self.assertFalse(result['all_stopping_checks_passed']);self.assertFalse(result['strict_state_comparisons_passed'])
    def test_missing_reference(self):
        label,result=self.fixture();row=result['records'].pop(label);candidate=label[:-1]+'1';row['repeat']=1;result['records'][candidate]=row
        b.finalize(result,[candidate]);self.assertEqual(result['comparison_missing_reference'],[candidate]);self.assertFalse(result['strict_state_comparisons_passed'])
    def test_empty(self):
        _,result=self.fixture();result['records']={};b.finalize(result,[])
        self.assertFalse(result['all_stopping_checks_passed']);self.assertFalse(result['strict_state_comparisons_passed'])
    def test_artifact_changed(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'state.npz';path.write_bytes(b'original')
            result=dict(records=dict(worker=dict(state=str(path),state_sha256=b.digest(path))),failures={})
            b.verify_artifacts(result);path.write_bytes(b'changed')
            with self.assertRaises(RuntimeError):b.verify_artifacts(result)
    def test_loaded_runtime_mismatch(self):
        label,result=self.fixture();row=result['records'][label];row['runtime_images']={'unexpected.so':'bad'}
        c=b.protocol_checks(row,result['manifest']['variants'][0],result['images']['reference'],result['configurations']['40:80:16']);self.assertFalse(c['passed']);self.assertFalse(c['runtime_images'])
    def test_physical_config_mismatch(self):
        _,result=self.fixture();row=next(iter(result['records'].values()));row['config']['hole'][0]['spin'][0]+=.001
        c=b.protocol_checks(row,result['manifest']['variants'][0],result['images']['reference'],result['configurations']['40:80:16']);self.assertFalse(c['passed']);self.assertFalse(c['configuration'])

if __name__=='__main__':unittest.main()
