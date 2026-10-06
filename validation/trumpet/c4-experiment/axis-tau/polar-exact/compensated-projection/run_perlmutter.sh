#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/axis-tau-compensated-v1"
mkdir "$PUNCTURE_RUN"
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_ROOT/source-axis-tau-compensated" -B "$PUNCTURE_ROOT/build-axis-tau-compensated" \
  -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_REGULARITY_CAP=4 \
  -DHISPID_RADIAL_STRETCH=.2 -DHISPID_ANGULAR_STRETCH=2. -DHISPID_NEWTON_BACKTRACKS=10 -DHISPID_STABLE_SCALAR_SOURCE=OFF -DHISPID_MONOTONE_PRECONDITIONER=OFF -DHISPID_AXIS_TAU=ON -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_COMPENSATED_MODAL_PROJECTION=ON \
  -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" \
  -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc \
  -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-axis-tau-compensated" --target HiSpID test_modal_projection -j4 > "$PUNCTURE_RUN/build.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" "$PUNCTURE_ROOT/source-axis-tau-compensated/src/HiSpID_solver.cpp" > "$PUNCTURE_RUN/images-sources.sha256"
"$PUNCTURE_ROOT/build-axis-tau-compensated/test_modal_projection" > "$PUNCTURE_RUN/cancellation-control.txt"
export HISPID_PROBE_LINEAR_ACTION=1 HISPID_TRACE_NEWTON=1
python3 "$PUNCTURE_ROOT/source-axis-tau-compensated/validation/run_trumpet_pilot.py" \
 --library "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" \
 --case moderate --n 240 --npolar 480 --nphi 16 --krylov-restart 64 \
 --memory-mib 65536 --linear-rtol .1 \
 --initial "$PUNCTURE_ROOT/axis-tau-polar-sequence-v3/moderate240/diagnostic.checkpoint" \
 --output "$PUNCTURE_RUN/moderate240" > "$PUNCTURE_RUN/solve240.log" 2>&1
