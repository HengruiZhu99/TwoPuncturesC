#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE="$PUNCTURE_ROOT/source-modal-relaxation-20261018"
PUNCTURE_BUILD="$PUNCTURE_ROOT/build-modal-relaxation-20261018"
PUNCTURE_RUN="$PUNCTURE_ROOT/modal-relaxation-build-v2"
mkdir "$PUNCTURE_RUN"
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_SOURCE" -B "$PUNCTURE_BUILD" \
 -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_REGULARITY_CAP=4 \
 -DHISPID_AXIS_TAU=ON -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_COMPENSATED_MODAL_PROJECTION=ON \
 -DHISPID_OVERINTEGRATED_GALERKIN=OFF -DHISPID_DAMPED_MODAL_RELAXATION=ON \
 -DHISPID_BATCHED_MODAL_TRANSFER=ON -DPUNCTURES_KOKKOS=ON \
 -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" \
 -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc \
 -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_BUILD" --target HiSpID -j2 > "$PUNCTURE_RUN/build.log" 2>&1
sha256sum "$PUNCTURE_BUILD/libHiSpID.so" "$PUNCTURE_BUILD/libTwoPunctures.so" > "$PUNCTURE_RUN/images.sha256"
