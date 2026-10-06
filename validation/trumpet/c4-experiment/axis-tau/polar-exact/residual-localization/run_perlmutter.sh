#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/axis-tau-residual-localization-v1"
mkdir "$PUNCTURE_RUN"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
printf '%s  %s\n' 2cfec4d71e54542add9805ee544dc73fbe4f54c3af8948e8ba1aabfabc61cdf1 "$PUNCTURE_ROOT/build-axis-tau-lgmres/libHiSpID.so" | sha256sum -c -
python3 "$PUNCTURE_ROOT/source-axis-tau-lgmres/validation/inspect_trumpet_residual.py" \
 --library "$PUNCTURE_ROOT/build-axis-tau-lgmres/libHiSpID.so" \
 --source-library "$PUNCTURE_ROOT/build-axis-tau-polar/libHiSpID.so" \
 --checkpoint "$PUNCTURE_ROOT/axis-tau-polar-sequence-v2/moderate224/diagnostic.checkpoint" \
 --output "$PUNCTURE_RUN/residual224.json" > "$PUNCTURE_RUN/residual224.log" 2>&1
python3 "$PUNCTURE_ROOT/source-axis-tau-lgmres/validation/inspect_trumpet_residual.py" \
 --library "$PUNCTURE_ROOT/build-axis-tau-lgmres/libHiSpID.so" \
 --source-library "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" \
 --checkpoint "$PUNCTURE_ROOT/axis-tau-compensated-v3/moderate240/diagnostic.checkpoint" \
 --output "$PUNCTURE_RUN/residual240.json" > "$PUNCTURE_RUN/residual240.log" 2>&1
