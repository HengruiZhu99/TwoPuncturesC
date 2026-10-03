"""Worker dispatch must reach its isolated sampler process entry point."""
import sys
from pathlib import Path
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'python'), str(ROOT/'examples'), str(ROOT/'validation')]
import portable_sampler_migration


class PortableSamplerDispatchTests(unittest.TestCase):
    def test_worker_command_reaches_sampler_entry_point(self):
        arguments = ['portable_sampler_migration.py', '--checkpoint', 'data.checkpoint',
                     '--worker-library', 'producer.so', '--worker-sha', 'a'*64,
                     '--worker-output', 'producer.npz']
        with patch.object(sys, 'argv', arguments), patch.object(
                portable_sampler_migration, 'worker') as entry_point:
            self.assertEqual(portable_sampler_migration.main(), 0)
            entry_point.assert_called_once()
            parsed = entry_point.call_args.args[0]
            self.assertEqual(parsed.worker_library, 'producer.so')
            self.assertEqual(parsed.worker_output, 'producer.npz')


if __name__ == '__main__':
    unittest.main()
