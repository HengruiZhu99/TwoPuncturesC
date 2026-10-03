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


def expected_workers(result):
    if set(result.get('records',{}))&set(result.get('failures',{})):
        raise ValueError('retained success/failure worker identities overlap')
    if 'expected_worker_ids' in result:
        workers=result['expected_worker_ids']
        if not isinstance(workers,list) or not workers or len(set(workers))!=len(workers):
            raise ValueError('explicit worker inventory must be nonempty and unique')
        if (set(result['records'])|set(result['failures']))-set(workers):
            raise ValueError('unexpected retained worker identity')
        return len(workers)
    binding=result['binding']
    return len(binding['grids'])*binding['repeats']*len(binding['systems'])*len(binding['methods'])*len(result['manifest']['variants'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job-id', required=True)
    parser.add_argument('--results', required=True)
    parser.add_argument('--finalize-margin-seconds', type=int, default=900)
    parser.add_argument('--wait-for-epoch', action='store_true', help='allow up to180 seconds for a newly resumed coordinator to publish this allocation')
    args = parser.parse_args()
    if not args.job_id.isdecimal() or args.finalize_margin_seconds < 60:
        parser.error('numeric job ID and at least 60 seconds finalization margin required')
    path = Path(args.results).resolve(strict=True)
    raw = path.parent/'raw'/path.stem
    if not raw.is_dir():
        raise FileNotFoundError(raw)
    observed_epoch = False
    epoch_deadline = time.monotonic() + 180
    while True:
        result = json.loads(path.read_text())  # atomic coordinator checkpoints
        epoch = result['allocation_epochs'][-1]
        if epoch['environment']['SLURM_JOB_ID'] != args.job_id:
            if args.wait_for_epoch and not observed_epoch:
                if time.monotonic() > epoch_deadline:
                    raise RuntimeError('Resumed coordinator did not publish its allocation epoch')
                time.sleep(2)
                continue
            print('Allocation epoch changed; observer exits without action', flush=True)
            return
        observed_epoch = True
        binding = result['binding']
        expected = expected_workers(result)
        if len(result['records']) + len(result['failures']) == expected:
            print('All workers recorded; finalization remains with coordinator', flush=True)
            return
        listing = subprocess.check_output(['squeue', '-h', '-j', args.job_id, '-o', '%L'], text=True, timeout=15).strip()
        if not listing:
            raise RuntimeError('Allocation disappeared before a clean checkpoint')
        remaining = duration(listing)
        # A worker may start just after the last poll. Include both the30s
        # sleep and the bounded15s scheduler query in the safety threshold.
        threshold = binding['timeout'] + args.finalize_margin_seconds + 45
        if remaining <= threshold:
            (raw/'STOP_AFTER_WORKER').touch()
            print(f'Requested between-worker checkpoint: {remaining}s left; worker timeout {binding["timeout"]}s plus margin {args.finalize_margin_seconds}s', flush=True)
            return
        time.sleep(30)


if __name__ == '__main__':
    main()
