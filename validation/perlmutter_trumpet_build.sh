#!/usr/bin/env bash
# Build only the chosen CUDA solve and CPU sampling paths, on a compute node.
set -euo pipefail
: "${SLURM_JOB_ID:?run inside an allocation}"
PUNCTURE_ROOT=${1:?fresh isolated run root}
PUNCTURE_KOKKOS=${2:?pinned Kokkos checkout}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_ROOT/job-id.txt"
hostname > "$PUNCTURE_ROOT/build-node.txt"
nvidia-smi -L > "$PUNCTURE_ROOT/gpu.txt"
for PUNCTURE_STAGE in sampler cuda; do
  PUNCTURE_FLAGS=(-DPUNCTURES_KOKKOS=OFF -DCMAKE_CXX_COMPILER="$NVCC_WRAPPER_DEFAULT_COMPILER")
  if [[ "$PUNCTURE_STAGE" == cuda ]]; then
    PUNCTURE_FLAGS=(-DPUNCTURES_KOKKOS=ON -DPUNCTURES_KOKKOS_SOURCE="$PUNCTURE_KOKKOS" -DCMAKE_CXX_COMPILER="$PUNCTURE_KOKKOS/bin/nvcc_wrapper" -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON -DKokkos_ARCH_AMPERE80=ON)
  fi
  cmake -S "$PUNCTURE_ROOT/source" -B "$PUNCTURE_ROOT/build-$PUNCTURE_STAGE" -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc "${PUNCTURE_FLAGS[@]}" > "$PUNCTURE_ROOT/configure-$PUNCTURE_STAGE.log" 2>&1
  /usr/bin/time -v cmake --build "$PUNCTURE_ROOT/build-$PUNCTURE_STAGE" --target HiSpID -j4 > "$PUNCTURE_ROOT/build-$PUNCTURE_STAGE.log" 2>&1
  sha256sum "$PUNCTURE_ROOT/build-$PUNCTURE_STAGE/libHiSpID.so" "$PUNCTURE_ROOT/build-$PUNCTURE_STAGE/libTwoPunctures.so" > "$PUNCTURE_ROOT/images-$PUNCTURE_STAGE.sha256"
done
printf 'completed\n' > "$PUNCTURE_ROOT/build-completed.txt"
