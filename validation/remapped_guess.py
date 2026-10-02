"""Explicit changed-map initial guesses; this module never accepts a solution.

Prepare in a process loading only the target library, then pass the JSON to
run_validation.py --initial-guess. The original record/library/raw vector
remain separately bound; every target grid still solves and is checked anew.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from hispid import Backend
from checkpoints import ROOT,PARAMETERIZATION,select_record,restore_payload
from configs import as_dict
from prolong import remap_modal,validate_modal_vector


def file_sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as source:
        for chunk in iter(lambda:source.read(1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()


def map_pair(maps):
    if not isinstance(maps,dict) or set(maps)!={'radial_stretch','angular_stretch'}:
        raise ValueError('explicit radial_stretch and angular_stretch required')
    radial,angular=maps['radial_stretch'],maps['angular_stretch']
    if not np.isfinite([radial,angular]).all() or not .001<=radial<=1 or not .1<=angular<=6:
        raise ValueError('unsupported collocation maps')
    return radial,angular


def modal_family(token):
    return token=='modal_P_C2prolate_mapped_v2' or token.startswith('modal_P_C2prolate_map_v3_')

def check_map_id(token,maps):
    actual=map_pair(maps)
    if token=='modal_P_C2prolate_mapped_v2':expected=(.2,2.)
    elif token.startswith('modal_P_C2prolate_map_v3_r'):
        try:expected=tuple(map(float,token.removeprefix('modal_P_C2prolate_map_v3_r').split('_k')))
        except ValueError:raise ValueError('malformed source map identifier') from None
    else:raise ValueError('supported mapped modal P identifier required')
    if expected!=actual:raise ValueError('source map identifier/parameters mismatch')


def load_guess(path,backend):
    """Verify both provenances and the target vector before any solve context."""
    payload=json.loads(Path(path).read_text())
    if (payload.get('kind')!='remapped_modal_initial_guess_v1' or payload.get('acceptance') is not False
            or payload.get('fresh_solve_required') is not True):
        raise ValueError('unaccepted initial-guess payload required')
    source=payload['source'];target=payload['target'];record=source['record']
    if (file_sha(source['library'])!=record['library_sha256']
            or file_sha(source['raw'])!=source['raw_sha256']):
        raise ValueError('source library/raw provenance mismatch')
    if (record['unknown_parameterization']!=PARAMETERIZATION
            or not modal_family(record['unknown_parameterization_id'])):
        raise ValueError('source continuous modal P family mismatch')
    if record['config']['n']!=record['resolution']:raise ValueError('source config/grid mismatch')
    check_map_id(record['unknown_parameterization_id'],record['collocation_maps'])
    with np.load(source['raw'],allow_pickle=False) as saved:
        validate_modal_vector(saved['unknowns'],record['resolution'])
    if (target['library_sha256']!=backend.library_sha256()
            or target['unknown_parameterization_id']!=backend.parameterization()
            or target['unknown_parameterization']!=backend.parameterization_description()
            or target['collocation_maps']!=backend.parameterization_maps()
            or target['unknown_parameterization']!=record['unknown_parameterization']):
        raise ValueError('target library/basis/maps mismatch')
    if not modal_family(target['unknown_parameterization_id']):raise ValueError('target modal P family required')
    raw=Path(path).resolve().parent/payload['vector_file']
    if file_sha(raw)!=payload['vector_sha256']:raise ValueError('initial-guess vector SHA mismatch')
    with np.load(raw,allow_pickle=False) as saved:values=saved['unknowns'].copy()
    map_pair(target['collocation_maps']);shape=target['config']['n']
    values=validate_modal_vector(values,shape)
    return payload,values


def validate_config(payload,config):
    from run_validation import SOLVER_CONTROLS
    actual=as_dict(config);target=payload['target']['config'];source=payload['source']['record']['config']
    if actual!=target:raise ValueError('first solve configuration must equal the declared initial-guess target')
    if any(actual[k]!=source[k] for k in actual if k not in SOLVER_CONTROLS):
        raise ValueError('changed-map initial guess must preserve all physical free data and frame inputs')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-library',required=True);parser.add_argument('--target-library',required=True)
    parser.add_argument('--case',required=True);parser.add_argument('--resolution',type=int);parser.add_argument('--nphi',type=int)
    parser.add_argument('--target-shape',required=True,help='radial:polar:Fourier')
    parser.add_argument('--krylov-restart',type=int,default=64);parser.add_argument('--memory-mib',type=int,default=8192)
    parser.add_argument('--output',required=True)
    args=parser.parse_args();out=Path(args.output).resolve();vector=out.with_suffix('.npz')
    if out==vector:raise ValueError('metadata output and npz vector must have distinct paths; use a .json output')
    if out.exists() or vector.exists():raise ValueError('initial-guess outputs already exist')
    record=select_record(args.case,args.resolution,args.nphi);source_library=Path(args.source_library).resolve(strict=True)
    if file_sha(source_library)!=record['library_sha256']:raise ValueError('source library SHA mismatch')
    # Loading the source image here would risk dyld reusing it for the target.
    backend=Backend(args.target_library)
    if (record['unknown_parameterization']!=PARAMETERIZATION or not modal_family(record['unknown_parameterization_id'])
            or backend.parameterization_description()!=PARAMETERIZATION or not modal_family(backend.parameterization())):
        raise ValueError('same continuous modal P family required')
    check_map_id(record['unknown_parameterization_id'],record['collocation_maps'])
    cfg,old=restore_payload(record,backend.config());cfg.n[:]=list(map(int,args.target_shape.split(':')))
    cfg.krylov_restart=args.krylov_restart;cfg.memory_limit_mib=args.memory_mib
    values=remap_modal(old,record['resolution'],list(cfg.n),map_pair(record['collocation_maps']),map_pair(backend.parameterization_maps()))
    source_raw=ROOT/f"validation/raw/{record['case']}_{record['resolution'][0]}_{record['resolution'][2]}.npz"
    payload=dict(kind='remapped_modal_initial_guess_v1',acceptance=False,fresh_solve_required=True,
        source=dict(library=str(source_library),raw=str(source_raw),raw_sha256=file_sha(source_raw),record=record),
        target=dict(library=str(backend.path),library_sha256=backend.library_sha256(),
            unknown_parameterization=backend.parameterization_description(),unknown_parameterization_id=backend.parameterization(),
            collocation_maps=backend.parameterization_maps(),config=as_dict(cfg)),vector_file=vector.name)
    validate_config(payload,cfg)
    out.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(vector,unknowns=values)
    payload['vector_sha256']=file_sha(vector);out.write_text(json.dumps(payload,indent=2)+'\n')
    checked,restored=load_guess(out,backend)
    if not np.array_equal(restored,values):raise ValueError('initial-guess serialization changed coefficients')
    print(json.dumps(dict(output=str(out),source_shape=record['resolution'],target_shape=list(cfg.n),
        source_maps=record['collocation_maps'],target_maps=backend.parameterization_maps(),
        coefficients=len(values),acceptance=False,fresh_solve_required=True)))


if __name__=='__main__':main()
