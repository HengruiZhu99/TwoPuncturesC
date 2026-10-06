"""Replay one Newton step from a retained checkpoint, without physical acceptance."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'python')]
from hispid import Backend
from checkpoint_export import read_checkpoint


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library', type=Path, required=True)
    p.add_argument('--source-library', type=Path, required=True)
    p.add_argument('--checkpoint', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--linear-rtol', type=float, default=.001)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    c, values, meta = read_checkpoint(a.checkpoint)
    if hashlib.sha256(a.source_library.read_bytes()).hexdigest() != meta['source_library_sha256']:
        raise ValueError('checkpoint source image mismatch')
    b = Backend(str(a.library.resolve()))
    if b.parameterization() != meta['parameterization']:
        raise ValueError('one-step replay requires the identical basis and maps')
    c.max_newton = 1
    os.environ['HISPID_TRACE_NEWTON'] = '1'
    result = dict(checkpoint=meta, library_sha256=b.loaded_sha256,
                  driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  linear_rtol=a.linear_rtol, binary_acceptance=False)
    with b.create(c, execution='kokkos', geometry='host') as s:
        s.set_unknowns(values)
        result['diagnostics'] = s.solve(krylov='gmres', linear_rtol=a.linear_rtol)
        result['linear_history'] = s.linear_history()
        result['work_statistics'] = s.work_statistics()
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
