"""Request a between-worker checkpoint before this Slurm allocation expires.

This administrative observer performs no numerical work and never terminates
the worker or allocation. Launch outside the immutable benchmark source tree.
"""
import argparse
import json
import subprocess
import time
from pathlib import Path


def duration(text):
    days, clock = text.split('-', 1) if '-' in text else ('0', text)
    parts = list(map(int, clock.split(':')))
    if len(parts) == 2:
        parts.insert(0, 0)
    if len(parts) != 3:
        raise ValueError('unknown Slurm time-left format: ' + text)
    return int(days)*86400 + parts[0]*3600 + parts[1]*60 + parts[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job-id', required=True)
    parser.add_argument('--results', required=True)
    parser.add_argument('--finalize-margin-seconds', type=int, default=900)
    args = parser.parse_args()
    if not args.job_id.isdecimal() or args.finalize_margin_seconds < 60:
        parser.error('numeric job ID and at least 60 seconds finalization margin required')
    path = Path(args.results).resolve(strict=True)
    raw = path.parent/'raw'/path.stem
    if not raw.is_dir():
        raise FileNotFoundError(raw)
    while True:
        result = json.loads(path.read_text())  # atomic coordinator checkpoints
        epoch = result['allocation_epochs'][-1]
        if epoch['environment']['SLURM_JOB_ID'] != args.job_id:
            print('Allocation epoch changed; observer exits without action', flush=True)
            return
        binding = result['binding']
        expected = len(binding['grids'])*binding['repeats']*len(binding['systems'])*len(binding['methods'])*len(result['manifest']['variants'])
        if len(result['records']) + len(result['failures']) == expected:
            print('All workers recorded; finalization remains with coordinator', flush=True)
            return
        listing = subprocess.check_output(['squeue', '-h', '-j', args.job_id, '-o', '%L'], text=True, timeout=15).strip()
        if not listing:
            raise RuntimeError('Allocation disappeared before a clean checkpoint')
        remaining = duration(listing)
        threshold = binding['timeout'] + args.finalize_margin_seconds
        if remaining <= threshold:
            (raw/'STOP_AFTER_WORKER').touch()
            print(f'Requested between-worker checkpoint: {remaining}s left; worker timeout {binding["timeout"]}s plus margin {args.finalize_margin_seconds}s', flush=True)
            return
        time.sleep(30)


if __name__ == '__main__':
    main()
