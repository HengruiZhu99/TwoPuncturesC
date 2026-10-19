#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-low-single-final-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/lower-resolution-20261018-v1
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export HISPID_TEST_LIBRARY=$PUNCTURE_ROOT/build-focus-tau-sampler/libHiSpID.so
python3 "$PUNCTURE_SOURCE/tests/test_single_hole.py" > "$PUNCTURE_RUN/single-tests.log" 2>&1
python3 "$PUNCTURE_SOURCE/validation/single_hole_study.py" \
 --library "$HISPID_TEST_LIBRARY" --output "$PUNCTURE_RUN/single-hole" > "$PUNCTURE_RUN/single-hole.log" 2>&1
python3 "$PUNCTURE_SOURCE/validation/compare_single_native.py" \
 --library "$PUNCTURE_ROOT/build-by-capacity-20261018/libTwoPunctures.so" \
 --output "$PUNCTURE_RUN/single-native.json" > "$PUNCTURE_RUN/single-native.log" 2>&1
for PUNCTURE_CASE in moderate paper-binary; do
 python3 "$PUNCTURE_SOURCE/validation/binary_reference_norm.py" \
  --study "$PUNCTURE_RUN/$PUNCTURE_CASE" \
  --by-library "$PUNCTURE_ROOT/build-by-capacity-20261018/libTwoPunctures.so" \
  --hispid-library "$PUNCTURE_ROOT/build-focus-dd-cuda/libHiSpID.so" \
  --output "$PUNCTURE_RUN/$PUNCTURE_CASE-reference-norm.json" > "$PUNCTURE_RUN/$PUNCTURE_CASE-reference-norm.log" 2>&1
done
# The first moderate BY grid's bulk observer already sits at its FD6 floor. This
# necessary coarse prefix resolves the truncation regime without repeating
# any of the completed32/48/64 solves or changing their equations/observers.
python3 "$PUNCTURE_ROOT/source-final-20261018/validation/lower_resolution_study.py" \
 --hispid-library "$PUNCTURE_ROOT/build-focus-dd-cuda/libHiSpID.so" \
 --by-library "$PUNCTURE_ROOT/build-by-capacity-20261018/libTwoPunctures.so" \
 --case moderate --levels 8,12,16 --threads 16 \
 --reference-hispid "$PUNCTURE_ROOT/focus-fd6-v1/hispid224/physical.npz" \
 --reference-by "$PUNCTURE_ROOT/focus-fd6-v1/by224/physical.npz" \
 --output "$PUNCTURE_RUN/moderate-prefix" > "$PUNCTURE_RUN/moderate-prefix.log" 2>&1
python3 "$PUNCTURE_SOURCE/validation/binary_reference_norm.py" \
 --study "$PUNCTURE_RUN/moderate-prefix" \
 --by-library "$PUNCTURE_ROOT/build-by-capacity-20261018/libTwoPunctures.so" \
 --hispid-library "$PUNCTURE_ROOT/build-focus-dd-cuda/libHiSpID.so" \
 --output "$PUNCTURE_RUN/moderate-prefix-reference-norm.json"
