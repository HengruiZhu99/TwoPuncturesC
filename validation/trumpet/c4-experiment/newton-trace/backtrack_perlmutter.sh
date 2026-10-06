#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/c4-backtrack-v1"
mkdir "$PUNCTURE_RUN"
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_ROOT/source-c4-backtrack" -B "$PUNCTURE_ROOT/build-c4-backtrack" \
  -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_REGULARITY_CAP=6 \
  -DHISPID_RADIAL_STRETCH=.2 -DHISPID_ANGULAR_STRETCH=2. -DHISPID_NEWTON_BACKTRACKS=24 \
  -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" \
  -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc \
  -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-c4-backtrack" --target HiSpID -j4 > "$PUNCTURE_RUN/build.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-c4-backtrack/libHiSpID.so" "$PUNCTURE_ROOT/source-c4-backtrack/src/HiSpID_solver.cpp" > "$PUNCTURE_RUN/images-sources.sha256"
export HISPID_TRACE_NEWTON=1
python3 "$PUNCTURE_ROOT/source-c4-backtrack/validation/run_trumpet_pilot.py" \
  --library "$PUNCTURE_ROOT/build-c4-backtrack/libHiSpID.so" \
  --case moderate --n 192 --npolar 384 --nphi 16 \
  --krylov-restart 64 --memory-mib 65536 --linear-rtol .001 \
  --initial "$PUNCTURE_ROOT/c4-experiment/moderate192/diagnostic.checkpoint" \
  --initial-source-library "$PUNCTURE_ROOT/build-c4-experiment-cuda/libHiSpID.so" \
  --output "$PUNCTURE_RUN/solve" > "$PUNCTURE_RUN/solve.log" 2>&1
