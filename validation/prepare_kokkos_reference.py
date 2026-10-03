"""Extract the frozen reference and declare its capacity-only benchmark patch."""
import argparse,json,shutil,tarfile
from pathlib import Path
from benchmark_bowen_york import digest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',required=True);p.add_argument('--destination',required=True)
    args=p.parse_args();destination=Path(args.destination).resolve()
    if destination.exists():raise FileExistsError(destination)
    destination.mkdir(parents=True)
    with tarfile.open(args.archive) as archive:archive.extractall(destination,filter='data')
    source=destination/'src/HiSpID_geometry.cpp';original=digest(source)
    text=source.read_text();old='c.memory_limit_mib>8192';new='c.memory_limit_mib>65536'
    if text.count(old)!=1:raise RuntimeError('frozen capacity guard changed')
    source.write_text(text.replace(old,new))
    shutil.copyfile(Path(__file__).resolve().parents[1]/'CMakeLists.txt',destination/'CMakeLists.txt')
    record=dict(reference_commit='25ca0649d2951537c1513fbcf4a7304698d84955',archive_sha256=digest(args.archive),
                capacity_only_patch=dict(file='src/HiSpID_geometry.cpp',before=old,after=new,original_sha256=original,patched_sha256=digest(source)),
                build_file_sha256=digest(destination/'CMakeLists.txt'),equations_and_stopping_changed=False,
                note='Existing aggregate allocation formula is unchanged. The benchmark explicitly requests32768MiB;128x256x28/restart64 was rejected under the original8192MiB cap.')
    (destination/'reference-build-provenance.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))
if __name__=='__main__':main()
