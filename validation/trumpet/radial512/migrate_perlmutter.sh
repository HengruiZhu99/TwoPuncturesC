#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/radial512/migration"
mkdir -p "$PUNCTURE_RUN"
export PYTHONPATH="$PUNCTURE_ROOT/source-radial512/python:$PUNCTURE_ROOT/workflow-radial512/validation"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=false
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
python3 "$PUNCTURE_ROOT/workflow-radial512/validation/portable_sampler_migration.py" --checkpoint "$PUNCTURE_ROOT/radial512/gamma10_512x320-v3/diagnostic.checkpoint" --producer-library "$PUNCTURE_ROOT/build-radial512-cuda/libHiSpID.so" --consumer-library "$PUNCTURE_ROOT/build-radial512-sampler/libHiSpID.so" --output "$PUNCTURE_RUN/proof.json" > "$PUNCTURE_RUN/migration.log" 2>&1
