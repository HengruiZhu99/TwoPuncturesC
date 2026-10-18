#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-by-capacity-20261018
PUNCTURE_BUILD=$PUNCTURE_ROOT/build-by-capacity-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/by-matched-capacity-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
cmake -S "$PUNCTURE_SOURCE" -B "$PUNCTURE_BUILD" -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper \
 -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" \
 > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_BUILD" -j8 --target TwoPunctures > "$PUNCTURE_RUN/build.log" 2>&1
python3 "$PUNCTURE_SOURCE/tests/test_by_capacity.py" --library "$PUNCTURE_BUILD/libTwoPunctures.so" --output "$PUNCTURE_RUN/budget-test.json"
sha256sum "$PUNCTURE_BUILD/libTwoPunctures.so" > "$PUNCTURE_RUN/images.sha256"
for PUNCTURE_N in 160 192 224; do
 python3 "$PUNCTURE_SOURCE/validation/compare_fixed_by.py" \
  --library "$PUNCTURE_BUILD/libTwoPunctures.so" \
  --hispid-result "$PUNCTURE_ROOT/qi-regularity-wide-v1/moderate$PUNCTURE_N/result.json" \
  --observer "$PUNCTURE_ROOT/qi-regularity-wide-v1/moderate$PUNCTURE_N/physical.npz" \
  --output "$PUNCTURE_RUN/by$PUNCTURE_N" > "$PUNCTURE_RUN/by$PUNCTURE_N.log" 2>&1
done
