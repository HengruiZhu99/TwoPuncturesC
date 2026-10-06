#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export PYTHONPATH="$PUNCTURE_ROOT/source-polar/python"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
PUNCTURE_RUN="$PUNCTURE_ROOT/extreme-pilots"
mkdir -p "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/spin99-job-id.txt"
python3 "$PUNCTURE_ROOT/workflow-polar/validation/run_trumpet_tolerance.py" --library "$PUNCTURE_ROOT/build-polar-cuda/libHiSpID.so" --case spin99 --n 160 --npolar 320 --nphi 32 --krylov-restart 64 --memory-mib 65536 --output "$PUNCTURE_RUN/spin99_160" > "$PUNCTURE_RUN/spin99.log" 2>&1
