"""Prerequisite mutation/compilation guards, without loading native images."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/name) for name in ('python','examples','validation')]
from run_extreme_kokkos import frozen_inputs,prerequisites,verify_hashes,source_floor_passed,verify_source_floor_arrays,completed_seed_controls
from check_far_source_floor import isolated_controls
from hispid import Config,Hole
from configs import as_dict


class ExtremePrerequisiteTests(unittest.TestCase):
    def test_seed_controls_reject_partial_or_duplicate_cases(self):
        seed=dict(completed=True,passed=True,cases=[dict(case=case,completed=True)
            for case in ('spin99','gamma10')])
        completed_seed_controls(seed)
        for partial in (seed|{'completed':False},seed|{'cases':seed['cases'][:1]},
            seed|{'cases':[seed['cases'][0],seed['cases'][0]]},
            seed|{'cases':[seed['cases'][0],seed['cases'][1]|{'completed':False}]}):
            with self.subTest(partial=partial),self.assertRaisesRegex(ValueError,'completed fresh'):
                completed_seed_controls(partial)

    def test_completed_failed_seed_controls_stay_failed(self):
        seed=dict(completed=True,passed=False,cases=[dict(case=case,completed=True,passed=False)
            for case in ('spin99','gamma10')])
        self.assertIs(completed_seed_controls(seed),seed)
        self.assertFalse(seed['passed'])
        self.assertTrue(all(not case['passed'] for case in seed['cases']))

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

    def source_floor(self):
        floor=dict(schema='hispid_source_floor_v3',extreme_controls=True,library_sha256='producer',
            residual_scaling='sin3_alpha_beta',execution='kokkos',compiled_execution='Cuda',seed_mass=.5,
            coordinate_separation=25.,weighted_far_source_bound=1e-14,far_radius_minimum=100.,
            exact_correction=0.,passed=True,records=[])
        for label,spin,velocity in isolated_controls(True):
            for active,sign in enumerate((1,-1)):
                cfg=Config();cfg.n[:]=[8,16,8];cfg.memory_limit_mib=32768
                cfg.hole[0]=Hole(0,(12.5,0,0));cfg.hole[1]=Hole(0,(-12.5,0,0))
                cfg.hole[active]=Hole(.5,(sign*12.5,0,0),spin*.25,-sign*velocity)
                floor['records'].append(dict(case=label,active_hole=active,resolution=list(cfg.n),config=as_dict(cfg),completed=True,
                    passed=True,weighted_far_source_linf=[1e-16]*4,physical_equivalent_far_linf=[1e-16]*4,
                    far_node_count=5,far_radius_range=[100.,1000.]))
        return floor

    def test_source_floor_requires_exact_target_and_complete_coverage(self):
        floor=self.source_floor();args=('producer',[[8,16,8]],25.)
        self.assertTrue(source_floor_passed(floor,*args))
        with self.assertRaises(ValueError):source_floor_passed(floor,'producer',[[16,32,8]],25.)
        with self.assertRaises(ValueError):source_floor_passed(floor,'producer',[[8,16,8]],26.)
        row=next(r for r in floor['records'] if r['case']=='spin99' and r['active_hole']==0)
        row['config']['hole'][0]['spin'][2]*=.95/.99
        with self.assertRaises(ValueError):source_floor_passed(floor,*args)

    def test_source_floor_requires_both_target_focus_velocity_signs(self):
        floor=self.source_floor();args=('producer',[[8,16,8]],25.)
        row=next(r for r in floor['records'] if r['case']=='gamma10' and r['active_hole']==0)
        row['config']['hole'][0]['velocity'][0]*=-1
        with self.assertRaises(ValueError):source_floor_passed(floor,*args)
        floor=self.source_floor();floor['records']=[r for r in floor['records'] if r['active_hole']==0]
        with self.assertRaises(ValueError):source_floor_passed(floor,*args)

    def test_numerical_floor_failure_remains_failed_despite_saved_pass_flag(self):
        floor=self.source_floor();args=('producer',[[8,16,8]],25.)
        floor['records'][0]['weighted_far_source_linf'][0]=2e-14
        self.assertFalse(source_floor_passed(floor,*args))
        floor=self.source_floor();floor['passed']=False
        self.assertFalse(source_floor_passed(floor,*args))
        floor=self.source_floor();floor['records'][0]['physical_equivalent_far_linf'][0]=float('nan')
        self.assertFalse(source_floor_passed(floor,*args))

    def test_source_floor_summaries_and_zero_state_are_recomputed(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'floor.npz';xyz=np.c_[np.arange(64)+100.,np.zeros((64,2))]
            def write(unknowns):
                np.savez_compressed(path,xyz=xyz,weighted=np.zeros((64,4)),physical_equivalent=np.zeros((64,4)),
                    unknowns=unknowns,far_mask=np.ones(64,dtype=bool))
                return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            row=dict(resolution=[4,4,4],far_node_count=64,far_radius_range=[100.,163.],
                weighted_far_source_linf=[0.]*4,physical_equivalent_far_linf=[0.]*4,raw_artifact=write(np.zeros(256)))
            floor=dict(records=[row],far_radius_minimum=100.)
            self.assertEqual(verify_source_floor_arrays(floor),{str(path):row['raw_artifact']['sha256']})
            row['far_node_count']=63
            with self.assertRaises(ValueError):verify_source_floor_arrays(floor)
            row['far_node_count']=64;row['raw_artifact']=write(np.ones(256))
            with self.assertRaisesRegex(ValueError,'zero-correction'):verify_source_floor_arrays(floor)


if __name__=='__main__':unittest.main()
