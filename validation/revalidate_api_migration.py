"""Fresh seed/moderate/covariance gates after a separately proved API migration.

Old checkpoints are read only through their bound old binary. The explicit
migration provides an initial guess, then the new binary solves and evaluates
all three moderate grids afresh. No old acceptance flags are relabeled.
"""
import argparse,copy,json,hashlib
from pathlib import Path
import numpy as np
from hispid import Backend,Config
from checkpoints import ROOT,restore_payload,select_record,library_sha
from run_validation import seeds,solve_case,REPORT,RAW
from check_covariance import run as covariance

p=argparse.ArgumentParser();p.add_argument('--old-library',required=True);p.add_argument('--new-library',required=True)
p.add_argument('--migration',default='validation/axis_api_migration.json')
p.add_argument('--source',default='moderate_far0_stable');p.add_argument('--label',default='moderate_sampler_api')
a=p.parse_args();new=Backend(a.new_library)
oldsha=hashlib.sha256(Path(a.old_library).resolve(strict=True).read_bytes()).hexdigest()
proof=json.loads((ROOT/a.migration).read_text())
if not proof['passed_off_axis_and_equations'] or not proof.get('isolated_processes') or not proof.get('loaded_images_verified') or proof['old_library_sha256']!=oldsha or proof['new_library_sha256']!=library_sha(new):
    raise ValueError('matching explicit equation/field migration proof is required')
archive=ROOT/'validation/before_sampler_api_results.json'
if not archive.exists():archive.write_text(REPORT.read_text())
source=select_record(a.source,24,12)
if source['library_sha256']!=oldsha:raise ValueError('initial guess source/build mismatch')
cfg,unknowns=restore_payload(source,new.config())
report=json.loads(REPORT.read_text());report['seeds']=seeds(new)
REPORT.write_text(json.dumps(report,indent=2)+'\n')
if not report['seeds']['passed']:raise SystemExit('new seed gate failed')
initial=copy.deepcopy(source);initial.update(case=a.label+'_migration_guess',library_sha256=library_sha(new),
                                           source_library_sha256=oldsha,migration_evidence=a.migration)
RAW.mkdir(exist_ok=True,parents=True)
np.savez_compressed(RAW/f"{initial['case']}_24_12.npz",unknowns=unknowns)
def factory(backend,n,nphi):
    result=Config.from_buffer_copy(cfg);result.n[:]=[n,n,nphi];return result
levels=[(24,12),(40,20),(56,28)]
result=solve_case(new,factory,levels,a.label,initial_record=initial)
result.update(stage='moderate',initial_guess_source_case=a.source,initial_guess_source_library_sha256=oldsha,migration_evidence=a.migration)
report=json.loads(REPORT.read_text());report[a.label]=result;REPORT.write_text(json.dumps(report,indent=2)+'\n')
if not result['passed']:raise SystemExit('new moderate gate failed')
report['covariance']=covariance(new,levels,a.label)
report=json.loads(REPORT.read_text())|{'covariance':report['covariance']}
REPORT.write_text(json.dumps(report,indent=2)+'\n')
raise SystemExit(0 if report['covariance']['passed'] else 1)
