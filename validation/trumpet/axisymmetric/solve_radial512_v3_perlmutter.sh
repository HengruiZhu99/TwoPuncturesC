#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/radial512"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source-radial512/python:$PUNCTURE_ROOT/workflow-radial512/validation"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/solve-v3-job.txt"
python3 "$PUNCTURE_ROOT/workflow-radial512/validation/run_trumpet_pilot.py" --library "$PUNCTURE_ROOT/build-radial512-cuda/libHiSpID.so" --case gamma10 --axisymmetric --n 512 --npolar 320 --nphi 8 --krylov-restart 64 --memory-mib 65536 --initial "$PUNCTURE_ROOT/axisymmetric/gamma10_256x320/diagnostic.checkpoint" --initial-source-library "$PUNCTURE_ROOT/build-axisym-cuda/libHiSpID.so" --output "$PUNCTURE_RUN/gamma10_512x320-v3" > "$PUNCTURE_RUN/solve-v3.log" 2>&1
