#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-final-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/lower-resolution-20261018-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python3 "$PUNCTURE_SOURCE/validation/lower_resolution_study.py" \
 --hispid-library "$PUNCTURE_ROOT/build-focus-dd-cuda/libHiSpID.so" \
 --by-library "$PUNCTURE_ROOT/build-by-capacity-20261018/libTwoPunctures.so" \
 --case moderate --levels 32,48,64 --threads 16 \
 --reference-hispid "$PUNCTURE_ROOT/focus-fd6-v1/hispid224/physical.npz" \
 --reference-by "$PUNCTURE_ROOT/focus-fd6-v1/by224/physical.npz" \
 --output "$PUNCTURE_RUN/moderate" > "$PUNCTURE_RUN/moderate.log" 2>&1
python3 "$PUNCTURE_SOURCE/validation/lower_resolution_study.py" \
 --hispid-library "$PUNCTURE_ROOT/build-focus-dd-cuda/libHiSpID.so" \
 --by-library "$PUNCTURE_ROOT/build-by-capacity-20261018/libTwoPunctures.so" \
 --case paper-binary --levels 16,24,32,48,72 --threads 16 \
 --output "$PUNCTURE_RUN/paper-binary" > "$PUNCTURE_RUN/paper-binary.log" 2>&1
