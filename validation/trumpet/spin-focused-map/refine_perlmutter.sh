#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source-polar/python"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_ROOT/spin-focused-map/refine-job-id.txt"
python3 "$PUNCTURE_ROOT/workflow-polar/validation/run_trumpet_tolerance.py" --library "$PUNCTURE_ROOT/build-spin-focused/libHiSpID.so" --case spin99 --n 192 --npolar 384 --nphi 32 --krylov-restart 64 --memory-mib 65536 --initial "$PUNCTURE_ROOT/spin-focused-map/spin99_160/diagnostic.checkpoint" --output "$PUNCTURE_ROOT/spin-focused-map/spin99_192" > "$PUNCTURE_ROOT/spin-focused-map/refine.log" 2>&1
