#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?run inside an allocation}"
PUNCTURE_ROOT=${1:?isolated run root}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_ROOT/consumer-job-id.txt"
# This single control checks device dispatch; the scientific suite is not repeated.
python3 "$PUNCTURE_ROOT/source/tests/check_trumpet_production.py" --library "$PUNCTURE_ROOT/build-cuda/libHiSpID.so" --output "$PUNCTURE_ROOT/device-control.json" > "$PUNCTURE_ROOT/device-control.log" 2>&1
cmake -S "$PUNCTURE_ROOT/athenak" -B "$PUNCTURE_ROOT/build-athenak" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ -DPROBLEM=z4c/hispid -DHISPID_ROOT="$PUNCTURE_ROOT/source" -DHISPID_LIBRARY_DIR="$PUNCTURE_ROOT/build-sampler" -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF -DKokkos_ENABLE_CUDA=OFF -DAthena_ENABLE_MPI=OFF -DAthena_ENABLE_OPENMP=OFF > "$PUNCTURE_ROOT/configure-athenak.log" 2>&1
/usr/bin/time -v cmake --build "$PUNCTURE_ROOT/build-athenak" -j8 > "$PUNCTURE_ROOT/build-athenak.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-athenak/src/athena" > "$PUNCTURE_ROOT/image-athenak.sha256"
printf 'completed\n' > "$PUNCTURE_ROOT/consumer-build-completed.txt"
