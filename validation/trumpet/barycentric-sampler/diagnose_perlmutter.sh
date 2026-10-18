#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE="$PUNCTURE_ROOT/source-barycentric-sampler-20261018"
PUNCTURE_BUILD="$PUNCTURE_ROOT/build-barycentric-sampler-20261018"
PUNCTURE_RUN="$PUNCTURE_ROOT/barycentric-sampler-diagnosis-v1"
mkdir "$PUNCTURE_RUN"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_SOURCE" -B "$PUNCTURE_BUILD" -DCMAKE_BUILD_TYPE=Release \
 -DHISPID_ROW_POWER=3 -DHISPID_REGULARITY_CAP=4 -DHISPID_AXIS_TAU=ON -DHISPID_EXACT_POLAR_TAU=ON \
 -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_BARYCENTRIC_SAMPLER=ON -DHISPID_GEOMETRIC_TRUMPET_INTERIOR=OFF \
 -DPUNCTURES_KOKKOS=OFF -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_BUILD" --target HiSpID test_trumpet_sampler -j4 > "$PUNCTURE_RUN/build.log" 2>&1
sha256sum "$PUNCTURE_BUILD/libHiSpID.so" > "$PUNCTURE_RUN/image.sha256"
"$PUNCTURE_BUILD/test_trumpet_sampler" > "$PUNCTURE_RUN/polynomial-control.txt"
python3 "$PUNCTURE_SOURCE/validation/trumpet/moderate-wide/diagnose_bulk.py" \
 --library "$PUNCTURE_BUILD/libHiSpID.so" --source-library "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" \
 --pilot "$PUNCTURE_ROOT/moderate-wide-sequence-v1/moderate224" --full-only \
 --output "$PUNCTURE_RUN/wide224.json" > "$PUNCTURE_RUN/wide224.log" 2>&1
python3 "$PUNCTURE_SOURCE/validation/trumpet/moderate-wide/diagnose_bulk.py" \
 --library "$PUNCTURE_BUILD/libHiSpID.so" --source-library "$PUNCTURE_ROOT/build-geometric-interior-20261018/libHiSpID.so" \
 --pilot "$PUNCTURE_ROOT/geometric-interior-sequence-v1/moderate192" --full-only \
 --output "$PUNCTURE_RUN/geometric192.json" > "$PUNCTURE_RUN/geometric192.log" 2>&1
