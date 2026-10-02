"""Build a frozen timing-only BY image; full-precision stdout changes no arithmetic."""
import argparse,json,shlex,subprocess
from pathlib import Path
from benchmark_bowen_york import ROOT,digest


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--objects',required=True);p.add_argument('--output-dir',required=True);args=p.parse_args()
    out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=False)
    source=out/'precision_newton.c';source.write_text((ROOT/'src/TP_Newton.c').read_text().replace('%10.3e','%.17e').replace('|F|=%e','|F|=%.17e'))
    wrapper=out/'timing.c';wrapper.write_text((ROOT/'validation/by_solver_timer.c').read_text().replace('"../src/TP_Newton.c"',json.dumps(str(source))))
    gsl_cflags=shlex.split(subprocess.check_output(['gsl-config','--cflags'],text=True));gsl_libs=shlex.split(subprocess.check_output(['gsl-config','--libs'],text=True))
    objects=sorted(p.resolve() for p in Path(args.objects).glob('*.o') if p.stem not in ('TP_Newton','TwoPuncturesRun'))
    if not objects or not (Path(args.objects)/'TP_Modal.o').exists():raise ValueError('native object directory missing modal solver')
    cmd=['gcc','-std=c99','-O3','-fPIC',*gsl_cflags,'-I'+str(ROOT/'include'),'-shared',str(wrapper),*map(str,objects),*gsl_libs,'-o',str(out/'libTwoPuncturesTimed.so')]
    subprocess.run(cmd,cwd=ROOT,check=True)
    record=dict(command=cmd,library_sha256=digest(out/'libTwoPuncturesTimed.so'),frozen_newton_sha256=digest(source),wrapper_sha256=digest(wrapper),source_sha256=digest(ROOT/'src/TP_Newton.c'),objects_sha256={str(p):digest(p) for p in objects},timing_only=True)
    (out/'build.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
if __name__=='__main__':main()
