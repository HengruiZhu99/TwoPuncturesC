#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/gamma-focused-map"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source-polar/python"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
python3 "$PUNCTURE_ROOT/workflow-polar/validation/run_trumpet_tolerance.py" --library "$PUNCTURE_ROOT/build-spin-focused/libHiSpID.so" --case gamma10 --n 160 --npolar 320 --nphi 8 --krylov-restart 64 --memory-mib 65536 --initial "$PUNCTURE_ROOT/extreme-pilots/gamma10_160/diagnostic.checkpoint" --initial-source-library "$PUNCTURE_ROOT/build-polar-cuda/libHiSpID.so" --output "$PUNCTURE_RUN/gamma10_160" > "$PUNCTURE_RUN/solve.log" 2>&1
