"""Bind performance variants to built native/runtime images and source hashes."""
import argparse,json,re,subprocess
from pathlib import Path
from benchmark_bowen_york import ROOT,digest

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True);p.add_argument('--output',required=True);args=p.parse_args();root=Path(args.root).resolve()
    variants=[];builds={}
    for build,space,thread_counts in [('reference','reference',[1]),('serial','Serial',[1]),('openmp','OpenMP',[1,2,4,8,16]),('cuda','Cuda',[16])]:
        path=root/('build-'+build);cache=path/'CMakeCache.txt'
        if not cache.is_file():raise FileNotFoundError(cache)
        dependencies={}
        for image in (path/'libHiSpID.so',path/'libTwoPunctures.so'):
            listing=subprocess.check_output(['ldd',str(image)],text=True)
            if 'not found' in listing:raise RuntimeError('unresolved dependency: '+listing)
            dependencies[str(image)]={str(Path(name).resolve(strict=True)):digest(Path(name).resolve(strict=True)) for name in re.findall(r'(?:=>\s+|^\s*)(/[^\s]+)\s+\(',listing,re.M)}
        runtime=sorted({f for files in dependencies.values() for f in files if 'libkokkos' in Path(f).name})
        if build!='reference' and not runtime:raise RuntimeError('shared Kokkos runtime missing')
        builds[build]=dict(cache_sha256=digest(cache),cache=cache.read_text(),images={str(f):digest(f) for f in (path/'libHiSpID.so',path/'libTwoPunctures.so')},resolved_dependencies=dependencies,runtime_images={f:digest(f) for f in runtime},tests=(root/('test-'+build+'.log')).read_text() if build!='reference' else None)
        for threads in thread_counts:
            label=build if build!='openmp' else f'openmp{threads}'
            variants.append(dict(id=label,execution='reference' if build=='reference' else 'kokkos',space=space,threads=threads,hispid_library=str(path/'libHiSpID.so'),by_library=str(path/'libTwoPunctures.so'),runtime_images=runtime,resolved_dependency_images={f:sha for files in dependencies.values() for f,sha in files.items()}))
    source={str(f.relative_to(ROOT)):digest(f) for directory in ('src','include','python','validation','tests') for f in (ROOT/directory).glob('*') if f.is_file() and f.suffix in ('.c','.cpp','.h','.hpp','.inc','.py','.sh')}
    source['CMakeLists.txt']=digest(ROOT/'CMakeLists.txt')
    attestation=json.loads((root/'kokkos-source-attestation.json').read_text())
    if attestation['git_head']!='6739bc623081648af9e752b616d9671527922cbf' or attestation['git_status_porcelain']:raise RuntimeError('Kokkos checkout is not the clean pinned commit')
    for filename,sha in attestation['tracked_files_sha256'].items():
        if digest(root/'kokkos'/filename)!=sha:raise RuntimeError('Kokkos attested file differs: '+filename)
    result=dict(variants=variants,builds=builds,source_sha256=source,reference_provenance=json.loads((root/'frozen-reference/reference-build-provenance.json').read_text()),kokkos_commit=attestation['git_head'],kokkos_source_attestation=attestation,note='Both CPU candidates and reference use GCC14 Release flags; CUDA uses the same host compiler through nvcc_wrapper. Native C remains C99. Complete caches and resolved images are recorded; no fast-math flag is requested.')
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print('Bound',len(variants),'variants to native/runtime/source hashes',flush=True)
if __name__=='__main__':main()
