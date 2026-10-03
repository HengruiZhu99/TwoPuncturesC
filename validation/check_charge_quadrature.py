"""Replay a bound checkpoint to check ADM quadrature and coordinate covariance."""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from hispid import Backend
from checkpoints import restore,restore_equivalent
from physical import charges, extrapolate
from native_loader import loaded_kokkos_images
from charge_checks import validate_refinements,radial_fits,angular_changes,qualify_refinement

p = argparse.ArgumentParser()
p.add_argument('--library', required=True)
p.add_argument('--case', default='moderate_far0_stable')
p.add_argument('--resolution', type=int, default=40)
p.add_argument('--nphi',type=int,help='required if the case has several angular grids at this radial resolution')
p.add_argument('--radii', default='200')
p.add_argument('--quadratures', default='12:24,20:40,32:64')
p.add_argument('--native-only', action='store_true', help='check native quadrature and radial fits without the independent/covariance callbacks')
p.add_argument('--compatibility-proof',help='explicit bitwise field/operator proof for read-only API migration; never accepts a binary')
p.add_argument('--output', help='save incremental JSON evidence relative to the worktree or at an absolute path')
p.add_argument('--results',help='separate study results JSON; default retains the historical validation/results.json')
p.add_argument('--raw-directory',help='separate study raw NPZ directory')
p.add_argument('--qualification-bound',type=float,help='explicit independent/angular/radial absolute bound; leaves whole binary acceptance false')
args = p.parse_args()
root = Path(__file__).resolve().parents[1]
source_path=Path(args.results) if args.results else root/'validation/results.json'
source_bytes=source_path.read_bytes();source_sha=hashlib.sha256(source_bytes).hexdigest()
data = json.loads(source_bytes)
candidates=[r for r in data[args.case]['records']
            if r['resolution'][0]==args.resolution and (args.nphi is None or r['resolution'][2]==args.nphi)]
if len(candidates)!=1:
    raise ValueError('select a unique checkpoint with --resolution and --nphi')
rec=candidates[0]
backend = Backend(args.library)
sha = backend.library_sha256()
proof_bytes=Path(args.compatibility_proof).read_bytes() if args.compatibility_proof else None
cfg,unknowns=restore_equivalent(backend,rec,json.loads(proof_bytes),args.raw_directory) if args.compatibility_proof else restore(backend,rec,args.raw_directory)
radii = [float(r) for r in args.radii.split(',')]
quadratures = [tuple(map(int, q.split(':'))) for q in args.quadratures.split(',')]
validate_refinements(radii,quadratures,args.qualification_bound)
if args.qualification_bound is not None:
    if not rec.get('raw_artifact_sha256'):raise ValueError('qualification requires retained raw artifact hashes')
raw_root=Path(args.raw_directory) if args.raw_directory else root/'validation/raw'
payload_path=raw_root/f"{rec['case']}_{rec['resolution'][0]}_{rec['resolution'][2]}.npz"
frozen_files={str(backend.path):sha}|backend.dependency_images|loaded_kokkos_images()
for original,value in rec.get('raw_artifact_sha256',{}).items():frozen_files[str(raw_root/Path(original).name)]=value
if not rec.get('raw_artifact_sha256'):frozen_files[str(payload_path)]=hashlib.sha256(payload_path.read_bytes()).hexdigest()
if proof_bytes is not None:frozen_files[str(Path(args.compatibility_proof))]=hashlib.sha256(proof_bytes).hexdigest()
axis = np.array([1., 2., 3.]); axis /= np.linalg.norm(axis)
W = np.array([[0, -axis[2], axis[1]],
              [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
Q = np.eye(3)+np.sin(.73)*W+(1-np.cos(.73))*(W@W)
offset = np.array([.7, -.2, .4])
def rotate_charge(q):
    return np.r_[q[0], Q@q[1:4], Q@q[4:7]]
output = dict(case=args.case, resolution=rec['resolution'],
              library_sha256=sha,unknown_parameterization_id=backend.parameterization(),
              checkpoint_source_library_sha256=rec['library_sha256'],compatibility_proof=args.compatibility_proof,
              collocation_maps=backend.parameterization_maps(),radii=radii,records=[],acceptance=False,
              results_source=str(source_path.resolve()),results_sha256=source_sha,
              solved_coordinate_covariance_verified=False,sampled_coordinate_transform_only=True)
output['bound_file_sha256']=frozen_files
if args.output and (root/args.output).exists():raise FileExistsError('use a fresh output to preserve quadrature evidence')
def save():
    if hashlib.sha256(source_path.read_bytes()).hexdigest()!=source_sha:
        raise ValueError('source study results changed during quadrature')
    if backend.library_sha256()!=sha or any(hashlib.sha256(Path(path).read_bytes()).hexdigest()!=value for path,value in frozen_files.items()):
        raise ValueError('bound quadrature image, proof or raw artifact changed')
    if args.output:
        path=root/args.output;path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(output,indent=2)+'\n')
save()
with backend.create_sampler(cfg) as solution:
    solution.set_unknowns(unknowns)  # Deliberately no solve().
    def rotated_sample(x):
        values = solution.sample((np.asarray(x)-offset)@Q)
        for key in ('gamma', 'Kij'):
            values[key] = np.einsum('ik,nkl,jl->nij', Q,
                values[key].reshape(-1, 3, 3), Q).reshape(-1, 9)
        return values
    previous = None
    for nt, np_ in quadratures:
        start=time.monotonic()
        native = np.array([solution.charges(r, ntheta=nt, nphi=np_)
                           for r in radii])
        if args.native_only:
            item = dict(ntheta=nt, nphi=np_, native_EPJ=native.tolist(),seconds=time.monotonic()-start)
            if previous is not None:
                item['quadrature_change_EPJ_linf'] = float(np.max(abs(native-previous)))
            previous = native
            if len(radii) >= 3:
                item['native_extrapolated_EPJ'] = extrapolate(radii,native).tolist()
            radial_fits(item,radii,native)
            angular_changes(item,output['records'][-1] if output['records'] else None)
            output['records'].append(item)
            save()
            print(json.dumps(item), flush=True)
            continue
        raw = np.array([charges(solution.sample, r, ntheta=nt, nphi=np_)
                        for r in radii])
        rotated = np.array([charges(rotated_sample, r, center=offset,
                                    ntheta=nt, nphi=np_) for r in radii])
        expected = np.array([rotate_charge(q) for q in raw])
        item = dict(ntheta=nt, nphi=np_, native_EPJ=native.tolist(),
                    independent_EPJ=raw.tolist(),
                    native_vs_independent_linf=float(np.max(abs(native-raw))),
                    centered_rotation_linf=float(np.max(abs(rotated-expected))))
        if previous is not None:
            item['quadrature_change_EPJ_linf'] = float(np.max(abs(raw-previous)))
        previous = raw
        if len(radii) >= 3:
            global_origin = np.array([charges(rotated_sample, r,
                ntheta=nt, nphi=np_) for r in radii])
            fit = extrapolate(radii, raw)
            expected_global = rotate_charge(fit)
            expected_global[4:7] += np.cross(offset, expected_global[1:4])
            item['independent_extrapolated_EPJ'] = fit.tolist()
            item['native_extrapolated_EPJ'] = extrapolate(radii,native).tolist()
            item['fixed_origin_extrapolation_linf'] = float(np.max(abs(
                extrapolate(radii, global_origin)-expected_global)))
        radial_fits(item,radii,raw)
        radial_fits(item,radii,native,'native_')
        if 'radial_fit_windows_EPJ' in item:
            item['native_vs_independent_radial_intercepts_linf']=float(np.max(abs(
                np.asarray(item['radial_fit_windows_EPJ'])-np.asarray(item['native_radial_fit_windows_EPJ']))))
        angular_changes(item,output['records'][-1] if output['records'] else None)
        output['records'].append(item)
        save()
        print(json.dumps(item), flush=True)
if args.qualification_bound is not None:
    output.update(qualify_refinement(output['records'],radii,args.qualification_bound))
    save()
print(json.dumps(output, indent=2))
