#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/axis-tau-polar-full-v1"
mkdir "$PUNCTURE_RUN"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export HISPID_TRACE_NEWTON=1
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
printf '%s  %s\n' 9489bca44c0bbd6ea20ae4996b4e2cdf225e84e3436ec80b253c4ffaf2b0c3e6 "$PUNCTURE_ROOT/build-axis-tau-polar/libHiSpID.so" | sha256sum -c -
python3 "$PUNCTURE_ROOT/source-axis-tau-polar/validation/run_trumpet_pilot.py" \
 --library "$PUNCTURE_ROOT/build-axis-tau-polar/libHiSpID.so" \
 --case moderate --n 192 --npolar 384 --nphi 16 \
 --krylov-restart 64 --memory-mib 65536 --linear-rtol .001 \
 --output "$PUNCTURE_RUN/moderate192" > "$PUNCTURE_RUN/solve.log" 2>&1
