#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-transition-focused-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/transition-focused-build-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
PUNCTURE_OPTIONS=(-DCMAKE_BUILD_TYPE=Release -DHISPID_INFINITY_POWER=1 -DHISPID_AXIS_TAU=ON -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_ROW_POWER=3 -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_BATCHED_MODAL_TRANSFER=ON -DHISPID_RADIAL_STRETCH=.1 -DHISPID_ANGULAR_STRETCH=3.)
cmake -S "$PUNCTURE_SOURCE" -B "$PUNCTURE_ROOT/build-transition-focused-cuda" \
 -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper \
 -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" \
 "${PUNCTURE_OPTIONS[@]}" > "$PUNCTURE_RUN/configure-cuda.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-transition-focused-cuda" -j8 --target HiSpID test_hispid_axis test_seed_regularity > "$PUNCTURE_RUN/build-cuda.log" 2>&1
"$PUNCTURE_ROOT/build-transition-focused-cuda/test_hispid_axis" > "$PUNCTURE_RUN/axis-test.log"
"$PUNCTURE_ROOT/build-transition-focused-cuda/test_seed_regularity" > "$PUNCTURE_RUN/seed-regularity-test.log"
cmake -S "$PUNCTURE_SOURCE" -B "$PUNCTURE_ROOT/build-transition-focused-sampler" \
 -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ "${PUNCTURE_OPTIONS[@]}" > "$PUNCTURE_RUN/configure-sampler.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-transition-focused-sampler" -j8 --target HiSpID > "$PUNCTURE_RUN/build-sampler.log" 2>&1
python3 "$PUNCTURE_SOURCE/tests/test_parallel_observer.py" --library "$PUNCTURE_ROOT/build-transition-focused-sampler/libHiSpID.so" --threads 16 --output "$PUNCTURE_RUN/observer-parity.json" > "$PUNCTURE_RUN/observer-parity.log" 2>&1
cmake -S "$PUNCTURE_ROOT/athenak-infinity2-20261018" -B "$PUNCTURE_ROOT/build-athenak-transition-focused" \
 -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ -DCMAKE_BUILD_TYPE=Release \
 -DPROBLEM=z4c/hispid -DAthena_ENABLE_OPENMP=ON -DAthena_ENABLE_MPI=OFF \
 -DHISPID_ROOT="$PUNCTURE_SOURCE" -DHISPID_LIBRARY_DIR="$PUNCTURE_ROOT/build-transition-focused-sampler" \
 > "$PUNCTURE_RUN/configure-athena.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-athenak-transition-focused" -j16 --target athena > "$PUNCTURE_RUN/build-athena.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-transition-focused-cuda/"lib*.so "$PUNCTURE_ROOT/build-transition-focused-sampler/"lib*.so "$PUNCTURE_ROOT/build-athenak-transition-focused/src/athena" > "$PUNCTURE_RUN/images.sha256"
