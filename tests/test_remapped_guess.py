"""Fail-closed provenance/configuration controls without numerical solves."""
import copy,json,tempfile,unittest
from pathlib import Path
import numpy as np
from unittest.mock import patch
from hispid import Config
from configs import as_dict
from checkpoints import PARAMETERIZATION
from remapped_guess import file_sha,load_guess,validate_config


class FakeBackend:
    def library_sha256(self):return 'a'*64
    def parameterization(self):return 'modal_P_C2prolate_map_v3_r0.2_k3'
    def parameterization_description(self):return PARAMETERIZATION
    def parameterization_maps(self):return dict(radial_stretch=.2,angular_stretch=3.)


class RemappedGuessTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.backend=FakeBackend();self.cfg=Config();self.cfg.n[:]=[8,8,8]
        self.cfg.hole[0].mass=.6;self.cfg.hole[1].mass=.4;self.cfg.hole[0].center[:]=[3,0,0];self.cfg.hole[1].center[:]=[-3,0,0]
        library=self.root/'source.so';library.write_bytes(b'preserved source library')
        self.source=self.root/'source.npz';self.vector=self.root/'guess.npz';self.metadata=self.root/'guess.json'
        self.values=np.arange(8*8*8*4,dtype=float)/100000
        np.savez_compressed(self.source,unknowns=self.values);np.savez_compressed(self.vector,unknowns=self.values)
        self.payload=dict(kind='remapped_modal_initial_guess_v1',acceptance=False,fresh_solve_required=True,
            source=dict(library=str(library),raw=str(self.source),raw_sha256=file_sha(self.source),record=dict(
                library_sha256=file_sha(library),unknown_parameterization=PARAMETERIZATION,
                unknown_parameterization_id='modal_P_C2prolate_mapped_v2',resolution=[8,8,8],
                collocation_maps=dict(radial_stretch=.2,angular_stretch=2.),config=as_dict(self.cfg))),
            target=dict(library_sha256=self.backend.library_sha256(),unknown_parameterization=PARAMETERIZATION,
                unknown_parameterization_id=self.backend.parameterization(),collocation_maps=self.backend.parameterization_maps(),config=as_dict(self.cfg)),
            vector_file=self.vector.name,vector_sha256=file_sha(self.vector))

    def save(self,payload=None):self.metadata.write_text(json.dumps(payload or self.payload))

    def test_valid_vector_roundtrip_and_free_data(self):
        self.save();payload,values=load_guess(self.metadata,self.backend)
        np.testing.assert_array_equal(values,self.values);validate_config(payload,self.cfg)

    def test_provenance_basis_map_acceptance_rejections(self):
        changes=[lambda p:p.update(acceptance=True),lambda p:p.update(fresh_solve_required=False),
            lambda p:p['source'].update(raw_sha256='0'*64),lambda p:p.update(vector_sha256='0'*64),
            lambda p:p['target'].update(library_sha256='0'*64),
            lambda p:p['target'].update(collocation_maps=dict(radial_stretch=.2,angular_stretch=2.)),
            lambda p:p['source']['record'].update(unknown_parameterization='legacy'),
            lambda p:p['source']['record'].update(collocation_maps=dict(radial_stretch=.2,angular_stretch=3.)),
            lambda p:p['source']['record'].update(resolution=[8,10,8])]
        for change in changes:
            with self.subTest(change=change):
                payload=copy.deepcopy(self.payload);change(payload);self.save(payload)
                with self.assertRaises(ValueError):load_guess(self.metadata,self.backend)

    def test_finite_vector_and_coefficient_count_rejections(self):
        for values in (self.values[:-1],np.full_like(self.values,np.nan)):
            for key,path in (('vector_sha256',self.vector),('raw_sha256',self.source)):
                with self.subTest(key=key,size=values.size):
                    np.savez_compressed(path,unknowns=values);payload=copy.deepcopy(self.payload)
                    (payload if key=='vector_sha256' else payload['source'])[key]=file_sha(path);self.save(payload)
                    with self.assertRaises(ValueError):load_guess(self.metadata,self.backend)
                    np.savez_compressed(path,unknowns=self.values)
                    # npz compression may change bytes across timestamps.
                    (self.payload if key=='vector_sha256' else self.payload['source'])[key]=file_sha(path)

    def test_changed_free_data_or_first_grid_rejected(self):
        payload=copy.deepcopy(self.payload);payload['source']['record']['config']['hole'][0]['center'][1]=.1
        with self.assertRaises(ValueError):validate_config(payload,self.cfg)
        payload=copy.deepcopy(self.payload);payload['target']['config']['n']=[10,10,8]
        with self.assertRaises(ValueError):validate_config(payload,self.cfg)

    def test_existing_case_or_raw_evidence_cannot_be_overwritten(self):
        import run_validation
        report=self.root/'results.json';report.write_text(json.dumps(dict(saved_case={})))
        with patch.object(run_validation,'REPORT',report),patch.object(run_validation,'RAW',self.root):
            with self.assertRaisesRegex(ValueError,'new label'):
                run_validation.solve_case(self.backend,None,[(8,8)],'saved_case',initial_guess=str(self.metadata))
            (self.root/'raw_case_8_8.npz').write_bytes(b'preserve earlier evidence')
            with self.assertRaisesRegex(ValueError,'new label'):
                run_validation.solve_case(self.backend,None,[(8,8)],'raw_case',initial_guess=str(self.metadata))

    def test_output_path_collision_rejected_before_loading(self):
        import remapped_guess
        args=['remapped_guess.py','--source-library','missing-source','--target-library','missing-target',
            '--case','unused','--target-shape','8:8:8','--output',str(self.root/'same.npz')]
        with patch('sys.argv',args),self.assertRaisesRegex(ValueError,'distinct paths'):remapped_guess.main()
        self.assertFalse((self.root/'same.npz').exists())

    def summary_records(self):
        records=[]
        for n,error in ((4,3e-8),(6,2e-8),(8,1e-8)):
            cfg=Config.from_buffer_copy(self.cfg);cfg.n[:]=[n,n,4]
            records.append(dict(case='summary',resolution=list(cfg.n),config=as_dict(cfg),
                library_sha256=self.backend.library_sha256(),unknown_parameterization=PARAMETERIZATION,
                unknown_parameterization_id=self.backend.parameterization(),collocation_maps=self.backend.parameterization_maps(),
                near=dict(H_rms=error,M_rms=error),bulk=dict(H_rms=error,M_rms=error),
                charges_extrapolated=[0.]*7,diagnostics=dict(status=0),min_metric_eigenvalue=1.,passed_local=True,passed_strict=True))
        np.savez_compressed(self.root/'summary_8_4.npz',unknowns=np.zeros(8*8*4*4))
        return records

    def summarize(self,records):
        import run_validation
        def factory(backend,n,p):
            cfg=Config.from_buffer_copy(self.cfg);cfg.n[:]=[n,n,p];return cfg
        with patch.object(run_validation,'RAW',self.root):
            return run_validation.solve_case(self.backend,factory,[],'summary',previous_records=records)

    def test_failed_earlier_solve_cannot_pass_strong_sequence(self):
        records=self.summary_records();self.assertTrue(self.summarize(records)['passed_strict'])
        records[0]['diagnostics']['status']=1;result=self.summarize(records)
        self.assertTrue(result['converges']);self.assertFalse(result['all_acceptance_solves_converged'])
        self.assertFalse(result['passed_strict'])

    def test_mixed_physical_sequence_and_repeated_grids_rejected(self):
        records=self.summary_records();records[0]['config']['hole'][0]['mass']=.7
        with self.assertRaisesRegex(ValueError,'physical free data'):self.summarize(records)
        records=self.summary_records()
        for r in records:r['resolution']=[8,8,4];r['config']['n']=[8,8,4]
        result=self.summarize(records);self.assertTrue(result['converges']);self.assertFalse(result['resolution_refines'])
        self.assertFalse(result['passed_strict'])

    def test_evaluate_only_cannot_relabel_another_stage(self):
        import run_validation
        report=self.root/'results.json';report.write_text(json.dumps(dict(mixed=dict(stage='moderate'))))
        args=['run_validation.py','--library','unused','--stage','spin95','--label','mixed','--resume','--evaluate-only']
        with patch('sys.argv',args),patch.object(run_validation,'REPORT',report),patch.object(run_validation,'Backend',return_value=self.backend):
            with self.assertRaisesRegex(ValueError,'another stage'):run_validation.main()


if __name__=='__main__':unittest.main()
