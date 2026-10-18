#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-focus2-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/focus2-consumer-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
PUNCTURE_OPTIONS=(-DCMAKE_BUILD_TYPE=Release -DHISPID_INFINITY_POWER=1 -DHISPID_AXIS_TAU=ON -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_ROW_POWER=3 -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_BATCHED_MODAL_TRANSFER=ON -DHISPID_FOCUS_TAU=ON -DHISPID_FOCUS_ORDER=2 -DHISPID_RADIAL_STRETCH=.2 -DHISPID_ANGULAR_STRETCH=2.)
cmake -S "$PUNCTURE_SOURCE" -B "$PUNCTURE_ROOT/build-focus2-sampler" \
 -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ "${PUNCTURE_OPTIONS[@]}" > "$PUNCTURE_RUN/configure-sampler.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-focus2-sampler" -j8 --target HiSpID > "$PUNCTURE_RUN/build-sampler.log" 2>&1
python3 "$PUNCTURE_SOURCE/tests/test_parallel_observer.py" --library "$PUNCTURE_ROOT/build-focus2-sampler/libHiSpID.so" --threads 16 --output "$PUNCTURE_RUN/observer-parity.json" > "$PUNCTURE_RUN/observer-parity.log" 2>&1
cmake -S "$PUNCTURE_ROOT/athenak-focus2-20261018" -B "$PUNCTURE_ROOT/build-athenak-focus2" \
 -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ -DCMAKE_BUILD_TYPE=Release \
 -DPROBLEM=z4c/hispid -DAthena_ENABLE_OPENMP=ON -DAthena_ENABLE_MPI=OFF \
 -DHISPID_ROOT="$PUNCTURE_SOURCE" -DHISPID_LIBRARY_DIR="$PUNCTURE_ROOT/build-focus2-sampler" \
 > "$PUNCTURE_RUN/configure-athena.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-athenak-focus2" -j16 --target athena > "$PUNCTURE_RUN/build-athena.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-focus2-sampler/"lib*.so "$PUNCTURE_ROOT/build-athenak-focus2/src/athena" > "$PUNCTURE_RUN/images.sha256"
