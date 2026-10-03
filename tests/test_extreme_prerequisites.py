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
from benchmark_kokkos import verify_artifacts


class ExtremePrerequisiteTests(unittest.TestCase):
    def complete_performance(self):
        grids=[[40,80,16],[80,160,16],[128,256,28]]
        variants=['reference','serial','openmp1','openmp2','openmp4','openmp8','openmp16','cuda']
        systems=['hispid','by'];methods=['gmres','bicgstab']
        records={f'{"_".join(map(str,grid))}_{variant}_{system}_{method}_{repeat}':
            dict(grid=grid,variant=variant,mode=system,krylov=method,repeat=repeat)
            for grid in grids for variant in variants for system in systems for method in methods for repeat in range(3)}
        failed=next(reversed(records));records.pop(failed)
        for label,row in records.items():
            for key in ('state','log','worker'):row.update({key:label+'.'+key,key+'_sha256':'0'*64})
        return dict(declared_performance_completed=True,expected_workers=288,completed_workers=288,
            binding=dict(grids=grids,systems=systems,methods=methods,repeats=3),
            manifest=dict(variants=[dict(id=v) for v in variants]),records=records,
            failures={failed:dict(log=failed+'.log',log_sha256='0'*64)},
            input_sha256={f'input_{"_".join(map(str,grid))}.json':'0'*64 for grid in grids})

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
        performance=self.complete_performance()
        receipt=dict(compilation_confirmed=True,compiler='mcp__codex_app__compile_latex_document',
                     performance_sha256=hashes['performance'],report_sha256=hashes['report'])
        prerequisites(performance,receipt,hashes)  # Retained failures do not become successful gates.
        for key,value in (('declared_performance_completed',False),('completed_workers',287),('expected_workers',287)):
            with self.subTest(key=key),self.assertRaises(ValueError):prerequisites(performance|{key:value},receipt,hashes)
        for key,value in (('compilation_confirmed',False),('compiler','unverified'),('performance_sha256','0'*64),('report_sha256','0'*64)):
            with self.subTest(key=key),self.assertRaises(ValueError):prerequisites(performance,receipt|{key:value},hashes)

    def test_completion_counts_cannot_hide_wrong_coverage_or_row_identity(self):
        hashes=dict(performance='performance',report='report')
        receipt=dict(compilation_confirmed=True,compiler='mcp__codex_app__compile_latex_document',
            performance_sha256='performance',report_sha256='report')
        performance=self.complete_performance();failed=next(iter(performance['failures']))
        performance['failures']={'wrong_key':{}}
        with self.assertRaises(ValueError):prerequisites(performance,receipt,hashes)
        performance=self.complete_performance();label=next(iter(performance['records']))
        performance['records'][label]['repeat']=3
        with self.assertRaisesRegex(ValueError,'record identity'):prerequisites(performance,receipt,hashes)
        performance=self.complete_performance();performance['manifest']['variants'][-1]={'id':'reference'}
        with self.assertRaises(ValueError):prerequisites(performance,receipt,hashes)
        performance=self.complete_performance();performance['records'][failed]=dict()
        with self.assertRaises(ValueError):prerequisites(performance,receipt,hashes)

    def test_retained_artifacts_rechecked_in_explicit_original_root(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);row={};paths=[]
            for key in ('state','log','worker'):
                path=root/key;path.write_bytes(key.encode());paths.append(path)
                row[key]=key;row[key+'_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
            input_path=root/'input';input_path.write_bytes(b'input');paths.append(input_path)
            result=dict(records={'worker':row},failures={},input_sha256={'input':hashlib.sha256(b'input').hexdigest()})
            verify_artifacts(result,root)
            for path in paths:
                before=path.read_bytes();path.write_bytes(before+b'changed')
                with self.subTest(artifact=path.name),self.assertRaises(RuntimeError):verify_artifacts(result,root)
                path.write_bytes(before)

    def test_missing_artifact_or_grid_input_witness_cannot_be_skipped(self):
        hashes=dict(performance='performance',report='report')
        receipt=dict(compilation_confirmed=True,compiler='mcp__codex_app__compile_latex_document',
            performance_sha256='performance',report_sha256='report')
        for key in ('state','state_sha256','log','log_sha256','worker','worker_sha256'):
            performance=self.complete_performance();next(iter(performance['records'].values())).pop(key)
            with self.subTest(missing=key),self.assertRaisesRegex(ValueError,'artifact witnesses'):
                prerequisites(performance,receipt,hashes)
        performance=self.complete_performance();next(iter(performance['failures'].values())).pop('log')
        with self.assertRaisesRegex(ValueError,'artifact witnesses'):prerequisites(performance,receipt,hashes)
        performance=self.complete_performance();performance['input_sha256']={}
        with self.assertRaisesRegex(ValueError,'input witnesses'):prerequisites(performance,receipt,hashes)
        performance=self.complete_performance();sha=performance['input_sha256'].pop('input_40_80_16.json')
        performance['input_sha256']['wrong_grid.json']=sha
        with self.assertRaisesRegex(ValueError,'input witnesses'):prerequisites(performance,receipt,hashes)

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

    def test_source_floor_binds_geometry_selection_and_native_witness(self):
        floor=self.source_floor();args=('producer',[[8,16,8]],25.)
        with self.assertRaises(ValueError):source_floor_passed(floor,*args,geometry='execution')
        floor['records'][0]['setup_statistics']=dict(geometry_execution=1,scalar_digits=53)
        with self.assertRaises(ValueError):source_floor_passed(floor,*args)
        floor['geometry']='execution'
        for row in floor['records']:row['setup_statistics']=dict(geometry_execution=1,scalar_digits=53)
        self.assertTrue(source_floor_passed(floor,*args,geometry='execution'))
        with self.assertRaises(ValueError):source_floor_passed(floor,*args)
        for witness in (None,dict(geometry_execution=0,scalar_digits=53),dict(geometry_execution=1,scalar_digits=64)):
            floor['records'][0]['setup_statistics']=witness
            with self.subTest(witness=witness),self.assertRaises(ValueError):source_floor_passed(floor,*args,geometry='execution')

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
