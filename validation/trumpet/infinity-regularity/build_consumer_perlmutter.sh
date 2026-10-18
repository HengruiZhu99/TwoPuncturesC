#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-infinity2-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/infinity2-consumer-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_SOURCE" -B "$PUNCTURE_ROOT/build-infinity2-sampler" \
 -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ -DCMAKE_BUILD_TYPE=Release \
 -DHISPID_INFINITY_POWER=2 -DHISPID_AXIS_TAU=ON -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_ROW_POWER=3 \
 -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_BATCHED_MODAL_TRANSFER=ON \
 > "$PUNCTURE_RUN/configure-native.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-infinity2-sampler" -j8 --target HiSpID > "$PUNCTURE_RUN/build-native.log" 2>&1
cmake -S "$PUNCTURE_ROOT/athenak-infinity2-20261018" -B "$PUNCTURE_ROOT/build-athenak-infinity2" \
 -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ -DCMAKE_BUILD_TYPE=Release \
 -DPROBLEM=z4c/hispid -DAthena_ENABLE_OPENMP=ON -DAthena_ENABLE_MPI=OFF \
 -DHISPID_ROOT="$PUNCTURE_SOURCE" -DHISPID_LIBRARY_DIR="$PUNCTURE_ROOT/build-infinity2-sampler" \
 > "$PUNCTURE_RUN/configure-athena.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-athenak-infinity2" -j16 --target athena > "$PUNCTURE_RUN/build-athena.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-infinity2-sampler/"lib*.so "$PUNCTURE_ROOT/build-athenak-infinity2/src/athena" > "$PUNCTURE_RUN/images.sha256"
