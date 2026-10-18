#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE="$PUNCTURE_ROOT/workflow-moderate-wide-20261018"
PUNCTURE_RUN="$PUNCTURE_ROOT/moderate-wide-pilot-v2"
mkdir "$PUNCTURE_RUN"
printf '%s  %s\n' 40d6a2165f07298427aa8acbca4d808748888eb8a579febe15815d02622f574c "$PUNCTURE_ROOT/build-axis-tau-compensated/libTwoPunctures.so" | sha256sum -c -
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores HISPID_TRACE_NEWTON=1
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
printf '%s  %s\n' 7398b4b254371c39409d511ebe0808fef5e8c45b9c1e6b1af2ccb52acecf31ea "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" | sha256sum -c -
python3 "$PUNCTURE_SOURCE/validation/run_trumpet_pilot.py" \
 --library "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" \
 --case moderate-wide --n 160 --npolar 320 --nphi 16 --krylov-restart 64 \
 --memory-mib 65536 --linear-rtol .1 \
 --output "$PUNCTURE_RUN/moderate160" > "$PUNCTURE_RUN/solve160.log" 2>&1
