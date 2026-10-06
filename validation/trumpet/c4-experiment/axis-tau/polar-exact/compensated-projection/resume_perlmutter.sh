#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/axis-tau-compensated-v2"
mkdir "$PUNCTURE_RUN"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
printf '%s  %s\n' 7398b4b254371c39409d511ebe0808fef5e8c45b9c1e6b1af2ccb52acecf31ea "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" | sha256sum -c -
printf '%s  %s\n' 9489bca44c0bbd6ea20ae4996b4e2cdf225e84e3436ec80b253c4ffaf2b0c3e6 "$PUNCTURE_ROOT/build-axis-tau-polar/libHiSpID.so" | sha256sum -c -
export HISPID_PROBE_LINEAR_ACTION=1 HISPID_TRACE_NEWTON=1
python3 "$PUNCTURE_ROOT/source-axis-tau-compensated/validation/run_trumpet_pilot.py" \
 --library "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" \
 --case moderate --n 240 --npolar 480 --nphi 16 --krylov-restart 64 \
 --memory-mib 65536 --linear-rtol .1 \
 --initial-source-library "$PUNCTURE_ROOT/build-axis-tau-polar/libHiSpID.so" \
 --initial "$PUNCTURE_ROOT/axis-tau-polar-sequence-v3/moderate240/diagnostic.checkpoint" \
 --output "$PUNCTURE_RUN/moderate240" > "$PUNCTURE_RUN/solve240.log" 2>&1
