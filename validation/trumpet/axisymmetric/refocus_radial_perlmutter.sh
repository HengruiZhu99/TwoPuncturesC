#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/gamma-radial-focused"
mkdir -p "$PUNCTURE_RUN"
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source-axisym/python"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_ROOT/source-axisym" -B "$PUNCTURE_ROOT/build-gamma-radial-focused" -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_RADIAL_STRETCH=.003 -DHISPID_ANGULAR_STRETCH=3.5 -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-gamma-radial-focused" --target HiSpID -j4 > "$PUNCTURE_RUN/build.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-gamma-radial-focused/libHiSpID.so" "$PUNCTURE_ROOT/source-axisym/src/HiSpID_solver.cpp" "$PUNCTURE_ROOT/source-axisym/src/HiSpID_symmetry.hpp" > "$PUNCTURE_RUN/images-sources.sha256"
python3 "$PUNCTURE_ROOT/workflow-axisym/validation/check_trumpet_axisymmetric.py" --library "$PUNCTURE_ROOT/build-gamma-radial-focused/libHiSpID.so" --execution kokkos --output "$PUNCTURE_RUN/control.json" > "$PUNCTURE_RUN/control.log" 2>&1
python3 "$PUNCTURE_ROOT/workflow-axisym/validation/run_trumpet_pilot.py" --library "$PUNCTURE_ROOT/build-gamma-radial-focused/libHiSpID.so" --case gamma10 --axisymmetric --n 256 --npolar 320 --nphi 8 --krylov-restart 64 --memory-mib 65536 --initial "$PUNCTURE_ROOT/axisymmetric/gamma10_256x320/diagnostic.checkpoint" --initial-source-library "$PUNCTURE_ROOT/build-axisym-cuda/libHiSpID.so" --output "$PUNCTURE_RUN/gamma10_256x320" > "$PUNCTURE_RUN/solve.log" 2>&1
