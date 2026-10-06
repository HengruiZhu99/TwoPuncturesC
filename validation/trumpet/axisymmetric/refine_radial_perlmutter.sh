#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source-axisym/python:$PUNCTURE_ROOT/workflow-axisym/validation"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_ROOT/axisymmetric/radial-job-id.txt"
python3 "$PUNCTURE_ROOT/workflow-axisym/validation/run_trumpet_pilot.py" --library "$PUNCTURE_ROOT/build-axisym-cuda/libHiSpID.so" --case gamma10 --axisymmetric --n 256 --npolar 320 --nphi 8 --krylov-restart 64 --memory-mib 65536 --initial "$PUNCTURE_ROOT/axisymmetric/gamma10_160/diagnostic.checkpoint" --output "$PUNCTURE_ROOT/axisymmetric/gamma10_256x320" > "$PUNCTURE_ROOT/axisymmetric/radial-refine.log" 2>&1
