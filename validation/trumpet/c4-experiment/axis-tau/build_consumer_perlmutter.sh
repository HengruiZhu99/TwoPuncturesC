#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/axis-tau-consumer-v1"
mkdir "$PUNCTURE_RUN"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_ROOT/source-axis-tau" -B "$PUNCTURE_ROOT/build-axis-tau-sampler" -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_RADIAL_STRETCH=.2 -DHISPID_ANGULAR_STRETCH=2. -DHISPID_REGULARITY_CAP=4 -DHISPID_AXIS_TAU=ON -DPUNCTURES_KOKKOS=OFF -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ > "$PUNCTURE_RUN/configure-sampler.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-axis-tau-sampler" --target HiSpID -j4 > "$PUNCTURE_RUN/build-sampler.log" 2>&1
cmake -S "$PUNCTURE_ROOT/athenak-affine" -B "$PUNCTURE_ROOT/build-athenak-axis-tau" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ -DPROBLEM=z4c/hispid -DHISPID_ROOT="$PUNCTURE_ROOT/source-axis-tau" -DHISPID_LIBRARY_DIR="$PUNCTURE_ROOT/build-axis-tau-sampler" -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=OFF -DAthena_ENABLE_MPI=OFF -DAthena_ENABLE_OPENMP=ON > "$PUNCTURE_RUN/configure-athenak.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-athenak-axis-tau" -j4 > "$PUNCTURE_RUN/build-athenak.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-axis-tau-sampler/libHiSpID.so" "$PUNCTURE_ROOT/build-axis-tau-sampler/libTwoPunctures.so" "$PUNCTURE_ROOT/build-athenak-axis-tau/src/athena" > "$PUNCTURE_RUN/images.sha256"
printf 'completed\n' > "$PUNCTURE_RUN/completed.txt"
