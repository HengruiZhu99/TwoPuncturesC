#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/spin-focused-map/consumer"
mkdir -p "$PUNCTURE_RUN"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_ROOT/source-polar" -B "$PUNCTURE_ROOT/build-focused-sampler" -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_RADIAL_STRETCH=.03 -DHISPID_ANGULAR_STRETCH=3.5 -DPUNCTURES_KOKKOS=OFF -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ > "$PUNCTURE_RUN/configure-sampler.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-focused-sampler" --target HiSpID -j8 > "$PUNCTURE_RUN/build-sampler.log" 2>&1
cmake -S "$PUNCTURE_ROOT/athenak-polar" -B "$PUNCTURE_ROOT/build-athenak-focused" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ -DPROBLEM=z4c/hispid -DHISPID_ROOT="$PUNCTURE_ROOT/source-polar" -DHISPID_LIBRARY_DIR="$PUNCTURE_ROOT/build-focused-sampler" -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=OFF -DAthena_ENABLE_MPI=OFF -DAthena_ENABLE_OPENMP=ON > "$PUNCTURE_RUN/configure-athenak.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-athenak-focused" -j8 > "$PUNCTURE_RUN/build-athenak.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-focused-sampler/libHiSpID.so" "$PUNCTURE_ROOT/build-focused-sampler/libTwoPunctures.so" "$PUNCTURE_ROOT/build-athenak-focused/src/athena" > "$PUNCTURE_RUN/images.sha256"
printf 'completed\n' > "$PUNCTURE_RUN/completed.txt"
