#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/axis-tau-bicgstab-probe-v1"
mkdir "$PUNCTURE_RUN"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
printf '%s  %s\n' 7398b4b254371c39409d511ebe0808fef5e8c45b9c1e6b1af2ccb52acecf31ea "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" | sha256sum -c -
python3 "$PUNCTURE_ROOT/source-axis-tau-compensated/validation/diagnose_trumpet_newton.py" \
 --library "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" \
 --source-library "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" \
 --checkpoint "$PUNCTURE_ROOT/axis-tau-compensated-v3/moderate240/diagnostic.checkpoint" \
 --krylov bicgstab --linear-rtol .1 --max-krylov 1200 --krylov-restart 64 \
 --output "$PUNCTURE_RUN/probe.json" > "$PUNCTURE_RUN/probe.log" 2>&1
