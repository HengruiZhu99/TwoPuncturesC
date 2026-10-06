#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source-polar/python:$PUNCTURE_ROOT/workflow-polar/validation"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_ROOT/spin-focused-map/replay-job-id.txt"
python3 "$PUNCTURE_ROOT/workflow-polar/validation/diagnose_trumpet_resolution.py" --library "$PUNCTURE_ROOT/build-spin-focused/libHiSpID.so" --checkpoint "$PUNCTURE_ROOT/spin-focused-map/spin99_192/diagnostic.checkpoint" --axes 0 1 --refined-shape 224 448 32 --output "$PUNCTURE_ROOT/spin-focused-map/directional-replay.json" > "$PUNCTURE_ROOT/spin-focused-map/directional-replay.log" 2>&1
