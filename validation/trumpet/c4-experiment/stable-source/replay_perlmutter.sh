#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/c4-stable-source-v2"
mkdir "$PUNCTURE_RUN"
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export HISPID_TRACE_NEWTON=1
python3 "$PUNCTURE_ROOT/source-c4-stable-source/validation/diagnose_trumpet_newton.py" \
  --library "$PUNCTURE_ROOT/build-c4-stable-source/libHiSpID.so" \
  --source-library "$PUNCTURE_ROOT/build-c4-experiment-cuda/libHiSpID.so" \
  --checkpoint "$PUNCTURE_ROOT/c4-experiment/moderate192/diagnostic.checkpoint" \
  --output "$PUNCTURE_RUN/newton.json" > "$PUNCTURE_RUN/newton.log" 2>&1
