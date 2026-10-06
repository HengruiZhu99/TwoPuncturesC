#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/c4-stable-source-v1"
mkdir "$PUNCTURE_RUN"
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_ROOT/source-c4-stable-source" -B "$PUNCTURE_ROOT/build-c4-stable-source" \
  -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_REGULARITY_CAP=6 \
  -DHISPID_RADIAL_STRETCH=.2 -DHISPID_ANGULAR_STRETCH=2. -DHISPID_NEWTON_BACKTRACKS=24 -DHISPID_STABLE_SCALAR_SOURCE=ON \
  -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" \
  -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc \
  -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-c4-stable-source" --target HiSpID -j4 > "$PUNCTURE_RUN/build.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-c4-stable-source/libHiSpID.so" "$PUNCTURE_ROOT/source-c4-stable-source/src/HiSpID_solver.cpp" > "$PUNCTURE_RUN/images-sources.sha256"
export HISPID_TRACE_NEWTON=1
python3 "$PUNCTURE_ROOT/source-c4-stable-source/validation/diagnose_trumpet_newton.py" \
  --library "$PUNCTURE_ROOT/build-c4-stable-source/libHiSpID.so" \
  --source-library "$PUNCTURE_ROOT/build-c4-experiment-cuda/libHiSpID.so" \
  --checkpoint "$PUNCTURE_ROOT/c4-experiment/moderate192/diagnostic.checkpoint" \
  --output "$PUNCTURE_RUN/newton.json" > "$PUNCTURE_RUN/newton.log" 2>&1
