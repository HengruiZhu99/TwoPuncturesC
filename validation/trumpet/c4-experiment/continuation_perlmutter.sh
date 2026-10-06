#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/c4-continuation-v1"
mkdir "$PUNCTURE_RUN"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export HISPID_TRACE_NEWTON=1
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
# Same unmodified preconditioner and extended-backtracking image as the last run.
python3 "$PUNCTURE_ROOT/source-c4-backtrack/validation/run_trumpet_pilot.py" \
  --library "$PUNCTURE_ROOT/build-c4-backtrack/libHiSpID.so" \
  --case moderate --n 12 --npolar 24 --nphi 16 \
  --krylov-restart 64 --memory-mib 65536 --linear-rtol .001 \
  --output "$PUNCTURE_RUN/coarse" > "$PUNCTURE_RUN/coarse.log" 2>&1
python3 - "$PUNCTURE_RUN/coarse/result.json" <<'PY'
import json,sys
r=json.load(open(sys.argv[1]))
if not r['completed'] or not r['diagnostics']['converged']:
    raise SystemExit('Coarse Newton solve did not converge; stopping continuation.')
PY
python3 "$PUNCTURE_ROOT/source-c4-backtrack/validation/run_trumpet_pilot.py" \
  --library "$PUNCTURE_ROOT/build-c4-backtrack/libHiSpID.so" \
  --case moderate --n 192 --npolar 384 --nphi 16 \
  --krylov-restart 64 --memory-mib 65536 --linear-rtol .001 \
  --initial "$PUNCTURE_RUN/coarse/diagnostic.checkpoint" \
  --output "$PUNCTURE_RUN/moderate192" > "$PUNCTURE_RUN/moderate192.log" 2>&1
