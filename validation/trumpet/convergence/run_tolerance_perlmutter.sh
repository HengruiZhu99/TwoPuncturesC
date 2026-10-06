#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export PYTHONPATH="$PUNCTURE_ROOT/source-polar/python"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_ROOT/convergence/tolerance-job-id.txt"
python3 "$PUNCTURE_ROOT/workflow-polar/validation/run_trumpet_tolerance.py" --library "$PUNCTURE_ROOT/build-polar-cuda/libHiSpID.so" --n 256 --npolar 512 --nphi 16 --krylov-restart 64 --memory-mib 65536 --tolerance 1e-14 --initial "$PUNCTURE_ROOT/convergence/moderate256/diagnostic.checkpoint" --output "$PUNCTURE_ROOT/convergence/moderate256tight" > "$PUNCTURE_ROOT/convergence/tolerance.log" 2>&1
