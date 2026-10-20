"""Workflow guards: native failure is terminal and saved progress is atomic."""
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'validation'))
from checkpointed_newton import solve_checkpointed


class Solution:
    config = SimpleNamespace(max_newton=1)

    def __init__(self, results, fail_write=False):
        self.results = iter(results)
        self.calls = 0
        self.fail_write = fail_write

    def solve(self, **options):
        assert options == dict(krylov='gmres', linear_rtol=.1)
        self.calls += 1
        return next(self.results)

    def linear_history(self):return [[0, .1, .09, 3]]
    def work_statistics(self):return dict(jvp_applications=4, preconditioner_applications=3)
    def write_checkpoint(self, path, acceptance):
        path.write_text(str(self.calls))
        if self.fail_write:raise OSError('interrupted write')
        return dict(path=str(path), acceptance=acceptance)


class CheckpointedNewton(unittest.TestCase):
    def run_case(self, results):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'progress.hispid';events=[];solution=Solution(results)
            result,attempts=solve_checkpointed(solution,path,8,'gmres',.1,lambda a,r:events.append(dict(r)))
            return solution.calls,result,attempts,events,path.read_text()

    def test_continues_only_native_newton_limit(self):
        calls,result,attempts,events,text=self.run_case([
            dict(converged=0,error='Newton iteration limit'),dict(converged=1,error='')])
        self.assertEqual(calls,2);self.assertEqual(len(attempts),2);self.assertEqual(text,'2')
        self.assertEqual([r['acceptance'] for r in events],['diagnostic','preliminary'])
        self.assertTrue(result['converged']);self.assertTrue(all(not r['path'].endswith('.pending') for r in events))

    def test_failed_linear_or_line_search_stops(self):
        for error in ('iteration limit or true residual above target','Newton line search failed'):
            calls,result,attempts,events,text=self.run_case([dict(converged=0,error=error)])
            self.assertEqual(calls,1);self.assertEqual(result['error'],error)
            self.assertEqual(events[0]['acceptance'],'diagnostic')

    def test_failed_write_preserves_previous_checkpoint(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'progress.hispid';path.write_text('previous')
            solution=Solution([dict(converged=0,error='Newton iteration limit')],fail_write=True)
            with self.assertRaises(OSError):solve_checkpointed(solution,path,8,'gmres',.1,lambda a,r:self.fail('published incomplete receipt'))
            self.assertEqual(path.read_text(),'previous')

if __name__=='__main__':unittest.main()
