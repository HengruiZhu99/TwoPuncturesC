"""Bind performance variants to built native/runtime images and source hashes."""
import argparse,json,re,subprocess
from pathlib import Path
from benchmark_bowen_york import ROOT,digest

def bind_build(path,tests=None):
    """Shared image/dependency binding for solve and supplemental setup runs."""
    path=Path(path).resolve(strict=True);cache=path/'CMakeCache.txt'
    if not cache.is_file():raise FileNotFoundError(cache)
    dependencies={}
    for image in (path/'libHiSpID.so',path/'libTwoPunctures.so'):
        listing=subprocess.check_output(['ldd',str(image)],text=True)
        if 'not found' in listing:raise RuntimeError('unresolved dependency: '+listing)
        dependencies[str(image)]={str(Path(name).resolve(strict=True)):digest(Path(name).resolve(strict=True)) for name in re.findall(r'(?:=>\s+|^\s*)(/[^\s]+)\s+\(',listing,re.M)}
    runtime=sorted({f for files in dependencies.values() for f in files if 'libkokkos' in Path(f).name})
    record=dict(cache_sha256=digest(cache),cache_path=str(cache),cache=cache.read_text(),images={str(f):digest(f) for f in (path/'libHiSpID.so',path/'libTwoPunctures.so')},resolved_dependencies=dependencies,runtime_images={f:digest(f) for f in runtime},tests=tests)
    return record,runtime,dependencies

def native_source(source):
    source=Path(source)
    files=[source/'CMakeLists.txt']+[f for directory in ('src','include','tests','validation') for f in (source/directory).glob('*') if f.is_file() and f.suffix in ('.c','.cpp','.h','.hpp','.inc')]
    return {str(f.relative_to(source)):digest(f) for f in files}

def setup_build_source(build,expected):
    """Bind imported new images to their actual CMake source and flags."""
    cache=build['cache']
    fields=dict(re.findall(r'^([^#/\r\n:=][^\r\n:=]*):[^\r\n=]+=(.*)$',cache,re.M))
    if fields.get('CMAKE_BUILD_TYPE')!='Release' or fields.get('PUNCTURES_KOKKOS')!='ON' or fields.get('PUNCTURES_BENCHMARK')!='ON' or fields.get('HISPID_ROW_POWER')!='3':
        raise ValueError('setup build must use Release/Kokkos/benchmark/cubic-row options')
    for key in ('CMAKE_C_FLAGS','CMAKE_CXX_FLAGS','CMAKE_C_FLAGS_RELEASE','CMAKE_CXX_FLAGS_RELEASE'):
        if re.search(r'(?:^|\s)(?:-ffast-math|-Ofast|--use_fast_math)(?:\s|$)',fields.get(key,'')):
            raise ValueError('setup build requests fast math')
    source=Path(fields['CMAKE_HOME_DIRECTORY']).resolve(strict=True)
    witness=native_source(source)
    if witness!=expected:raise ValueError('setup compiled native source differs from staged validation source')
    return dict(path=str(source),native_source_sha256=witness)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True);p.add_argument('--output',required=True)
    p.add_argument('--setup-build-map',help='JSON mapping serial/openmp/cuda to already compiled new build directories; emits a separate setup-only manifest')
    args=p.parse_args();root=Path(args.root).resolve()
    mapped=json.loads(Path(args.setup_build_map).read_text()) if args.setup_build_map else None
    if mapped is not None and set(mapped)!={'serial','openmp','cuda'}:raise ValueError('setup build map must name exactly serial/openmp/cuda')
    variants=[];builds={}
    for build,space,thread_counts in [('reference','reference',[1]),('serial','Serial',[1]),('openmp','OpenMP',[1,2,4,8,16]),('cuda','Cuda',[16])]:
        if mapped is not None and build=='reference':continue
        path=Path(mapped[build]).resolve(strict=True) if mapped is not None else root/('build-'+build)
        tests=None if mapped is not None or build=='reference' else (root/('test-'+build+'.log')).read_text()
        builds[build],runtime,dependencies=bind_build(path,tests)
        if build!='reference' and not runtime:raise RuntimeError('shared Kokkos runtime missing')
        if mapped is not None:
            builds[build]['compiled_source']=setup_build_source(builds[build],native_source(ROOT))
            builds[build]['actual_setup_tests_pending']=True
        for threads in thread_counts:
            label=build if build!='openmp' else f'openmp{threads}'
            variants.append(dict(id=label,execution='reference' if build=='reference' else 'kokkos',space=space,threads=threads,hispid_library=str(path/'libHiSpID.so'),by_library=str(path/'libTwoPunctures.so'),runtime_images=runtime,resolved_dependency_images={f:sha for files in dependencies.values() for f,sha in files.items()}))
    source={str(f.relative_to(ROOT)):digest(f) for directory in ('src','include','python','validation','tests') for f in (ROOT/directory).glob('*') if f.is_file() and f.suffix in ('.c','.cpp','.h','.hpp','.inc','.py','.sh')}
    source['CMakeLists.txt']=digest(ROOT/'CMakeLists.txt')
    attestation=json.loads((root/'kokkos-source-attestation.json').read_text())
    if attestation['git_head']!='6739bc623081648af9e752b616d9671527922cbf' or attestation['git_status_porcelain']:raise RuntimeError('Kokkos checkout is not the clean pinned commit')
    for filename,sha in attestation['tracked_files_sha256'].items():
        if digest(root/'kokkos'/filename)!=sha:raise RuntimeError('Kokkos attested file differs: '+filename)
    result=dict(variants=variants,builds=builds,source_sha256=source,kokkos_commit=attestation['git_head'],kokkos_source_attestation=attestation,note='CPU images use Release flags; CUDA uses the host compiler through nvcc_wrapper. Native C remains C99. Complete caches and resolved images are recorded; no fast-math flag is requested.')
    if mapped is None:result['reference_provenance']=json.loads((root/'frozen-reference/reference-build-provenance.json').read_text())
    else:result.update(scope='supplemental paired geometry setup',solve_campaign_replaced=False,setup_build_map_sha256=digest(args.setup_build_map),actual_setup_tests_pending=True)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print('Bound',len(variants),'variants to native/runtime/source hashes',flush=True)
if __name__=='__main__':main()
