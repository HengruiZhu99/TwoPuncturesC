#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/c4-experiment"
export PYTHONPATH="$PUNCTURE_ROOT/source-c4-experiment/python:$PUNCTURE_ROOT/source-c4-experiment/validation:$PUNCTURE_ROOT/source-c4-experiment/examples"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/linear-retry-job.txt"
# Diagnose the failed inexact Newton step; no physical free data or basis changes.
python3 "$PUNCTURE_ROOT/source-c4-experiment/validation/run_trumpet_pilot_linear.py" \
  --library "$PUNCTURE_ROOT/build-c4-experiment-cuda/libHiSpID.so" \
  --case moderate --n 192 --npolar 384 --nphi 16 \
  --krylov-restart 64 --memory-mib 65536 --linear-rtol .001 \
  --initial "$PUNCTURE_RUN/moderate192/diagnostic.checkpoint" \
  --output "$PUNCTURE_RUN/moderate192-linear001" > "$PUNCTURE_RUN/moderate192-linear001.log" 2>&1
