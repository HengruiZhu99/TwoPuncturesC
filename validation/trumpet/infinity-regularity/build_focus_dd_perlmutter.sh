#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-focus-dd-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/focus-dd-build-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_SOURCE" -B "$PUNCTURE_ROOT/build-focus-dd-cuda" -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper \
 -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" \
 -DHISPID_AXIS_TAU=ON -DHISPID_FOCUS_TAU=ON -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_ROW_POWER=3 \
 -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_BATCHED_MODAL_TRANSFER=ON \
 > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-focus-dd-cuda" -j8 --target HiSpID test_hispid_focus test_hispid_kokkos > "$PUNCTURE_RUN/build.log" 2>&1
"$PUNCTURE_ROOT/build-focus-dd-cuda/test_hispid_focus" > "$PUNCTURE_RUN/focus-test.log"
sha256sum "$PUNCTURE_ROOT/build-focus-dd-cuda/"lib*.so "$PUNCTURE_ROOT/build-focus-dd-cuda/test_hispid_kokkos" > "$PUNCTURE_RUN/images.sha256"
