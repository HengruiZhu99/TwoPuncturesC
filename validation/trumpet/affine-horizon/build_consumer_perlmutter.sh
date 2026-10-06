#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/affine-horizon/consumer"
mkdir -p "$PUNCTURE_RUN"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
/opt/cray/pe/gcc-native/14/bin/g++ -std=c++17 -O2 -I "$PUNCTURE_ROOT/athenak-affine/src" "$PUNCTURE_ROOT/athenak-affine/tst/test_suite/z4c/check_affine_horizon.cpp" -o "$PUNCTURE_RUN/check-affine"
"$PUNCTURE_RUN/check-affine" > "$PUNCTURE_RUN/algebra-control.log"
cmake -S "$PUNCTURE_ROOT/athenak-affine" -B "$PUNCTURE_ROOT/build-athenak-affine" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ -DPROBLEM=z4c/hispid -DHISPID_ROOT="$PUNCTURE_ROOT/source-polar" -DHISPID_LIBRARY_DIR="$PUNCTURE_ROOT/build-focused-sampler" -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=OFF -DAthena_ENABLE_MPI=OFF -DAthena_ENABLE_OPENMP=ON > "$PUNCTURE_RUN/configure-athenak.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-athenak-affine" -j4 > "$PUNCTURE_RUN/build-athenak.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-focused-sampler/libHiSpID.so" "$PUNCTURE_ROOT/build-focused-sampler/libTwoPunctures.so" "$PUNCTURE_ROOT/build-athenak-affine/src/athena" > "$PUNCTURE_RUN/images.sha256"
printf 'completed\n' > "$PUNCTURE_RUN/completed.txt"
