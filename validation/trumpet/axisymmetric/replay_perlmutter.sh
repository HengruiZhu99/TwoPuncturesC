#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source-axisym/python:$PUNCTURE_ROOT/workflow-axisym/validation"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_ROOT/axisymmetric/replay-job-id.txt"
python3 "$PUNCTURE_ROOT/workflow-axisym/validation/diagnose_trumpet_resolution.py" --library "$PUNCTURE_ROOT/build-axisym-cuda/libHiSpID.so" --checkpoint "$PUNCTURE_ROOT/axisymmetric/gamma10_160/diagnostic.checkpoint" --axes 0 1 --refined-shape 192 384 8 --output "$PUNCTURE_ROOT/axisymmetric/directional-replay.json" > "$PUNCTURE_ROOT/axisymmetric/directional-replay.log" 2>&1
