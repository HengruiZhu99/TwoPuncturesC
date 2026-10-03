"""Malformed interchange records must fail before any native image is loaded."""
import sys
import tempfile
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
from checkpoint_export import read_checkpoint,write_checkpoint
from hispid import Config,Hole


class PortableCheckpointTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.path=Path(self.directory.name)/'portable.txt'
        self.config=Config();self.config.n[:]=[4,4,4];self.config.memory_limit_mib=32768
        self.config.hole[0]=Hole(.5,(12.5,0,0),velocity=(-float(np.sqrt(.99)),0,0))
        self.config.hole[1]=Hole(.5,(-12.5,0,0),velocity=(float(np.sqrt(.99)),0,0))
        self.values=np.linspace(-.2,.3,256)
        write_checkpoint(self.path,self.config,self.values,'a'*64,'diagnostic')
        self.original=self.path.read_text()

    def tearDown(self):self.directory.cleanup()

    def test_typed_roundtrip_preserves_budget_fields_and_values(self):
        config,values,metadata=read_checkpoint(self.path)
        self.assertEqual(bytes(config),bytes(self.config))
        np.testing.assert_array_equal(values,self.values)
        self.assertEqual(metadata['acceptance'],'diagnostic')
        self.assertEqual(metadata['source_library_sha256'],'a'*64)

    def test_malformed_records_are_rejected(self):
        mutations=(self.original.replace('memory_limit_mib 32768','memory_limit_mib 2147483648'),
            self.original.replace('unknowns 256','unknowns 252'),
            self.original.replace('END\n',''),self.original+'unused\n',
            self.original.replace('hole0 0.5','hole0 nan'),
            self.original.replace('acceptance diagnostic','acceptance strong diagnostic'),
            self.original.replace('unknowns 256','memory_limit_mib 32768\nunknowns 256'),
            self.original.replace('modal_P_C2prolate_mapped_v2','modal_P_C2prolate_map_v3_rnan_k2'),
            self.original.replace('n 4 4 4','n 2147483647 2147483647 2147483647'))
        for text in mutations:
            with self.subTest(text=text[:70]):
                self.path.write_text(text)
                with self.assertRaises(ValueError):read_checkpoint(self.path)


if __name__=='__main__':unittest.main()
