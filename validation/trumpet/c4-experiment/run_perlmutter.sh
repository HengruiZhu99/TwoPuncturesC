#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/c4-experiment"
export PYTHONPATH="$PUNCTURE_ROOT/source-c4-experiment/python:$PUNCTURE_ROOT/source-c4-experiment/validation:$PUNCTURE_ROOT/source-c4-experiment/examples"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/production-job.txt"
python3 "$PUNCTURE_ROOT/source-c4-experiment/validation/check_regular_operators.py" --library "$PUNCTURE_ROOT/build-c4-experiment-cuda/libHiSpID.so" --levels 40 --nphi 16 --modes 5,6,8 --components 0,1 --regularity-cap 6 --execution kokkos --output "$PUNCTURE_RUN/operator-control.json" > "$PUNCTURE_RUN/operator-control.log" 2>&1
python3 "$PUNCTURE_ROOT/source-c4-experiment/validation/run_trumpet_pilot.py" --library "$PUNCTURE_ROOT/build-c4-experiment-cuda/libHiSpID.so" --case moderate --n 192 --npolar 384 --nphi 16 --krylov-restart 64 --memory-mib 65536 --output "$PUNCTURE_RUN/moderate192" > "$PUNCTURE_RUN/moderate192.log" 2>&1
