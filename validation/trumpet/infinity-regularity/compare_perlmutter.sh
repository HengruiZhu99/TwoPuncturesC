#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-infinity2-20261018
PUNCTURE_BUILD=$PUNCTURE_ROOT/build-infinity2-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/infinity2-by-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python3 "$PUNCTURE_SOURCE/validation/check_liu_seed.py" --library "$PUNCTURE_BUILD/libHiSpID.so" \
 --execution kokkos --output "$PUNCTURE_RUN/liu-seed-cuda.json" > "$PUNCTURE_RUN/liu.log" 2>&1
# One frozen-field check determines whether the remaining QI bulk floor
# comes from the held-fixed phi grid; there is no new nonlinear solve.
python3 "$PUNCTURE_SOURCE/validation/trumpet/moderate-wide/diagnose_phi_alias.py" \
 --library "$PUNCTURE_ROOT/build-stable-geometry-20261018/libHiSpID.so" \
 --checkpoint "$PUNCTURE_ROOT/stable-geometry-v1/moderate160/diagnostic.checkpoint" \
 --output "$PUNCTURE_RUN/qi-phi-alias.json" > "$PUNCTURE_RUN/phi-alias.log" 2>&1

# Reuse the baseline's immutable Cartesian observers. Neither the BY
# comparison nor an infinity-map trial can relocate these sample points.
for PUNCTURE_N in 160 192 224; do
 python3 "$PUNCTURE_SOURCE/validation/compare_fixed_by.py" \
  --library "$PUNCTURE_BUILD/libTwoPunctures.so" \
  --hispid-result "$PUNCTURE_ROOT/qi-regularity-wide-v1/moderate$PUNCTURE_N/result.json" \
  --observer "$PUNCTURE_ROOT/qi-regularity-wide-v1/moderate$PUNCTURE_N/physical.npz" \
  --output "$PUNCTURE_RUN/by$PUNCTURE_N" > "$PUNCTURE_RUN/by$PUNCTURE_N.log" 2>&1
done
