#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/radial512"
mkdir -p "$PUNCTURE_RUN"
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source-radial512/python"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_ROOT/source-radial512" -B "$PUNCTURE_ROOT/build-radial512-cuda" -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_RADIAL_STRETCH=.03 -DHISPID_ANGULAR_STRETCH=3.5 -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-radial512-cuda" --target HiSpID -j4 > "$PUNCTURE_RUN/build.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-radial512-cuda/libHiSpID.so" "$PUNCTURE_ROOT/source-radial512/src/HiSpID_solver.cpp" "$PUNCTURE_ROOT/source-radial512/src/HiSpID_symmetry.hpp" > "$PUNCTURE_RUN/images-sources.sha256"
export PYTHONPATH="$PUNCTURE_ROOT/source-radial512/python:$PUNCTURE_ROOT/workflow-axisym/validation"
python3 "$PUNCTURE_ROOT/check_regular_operators_radial512.py" --library "$PUNCTURE_ROOT/build-radial512-cuda/libHiSpID.so" --levels 512 --npolar 16 --nphi 8 --modes 0,1 --components 0,1 --memory-mib 4096 --output "$PUNCTURE_RUN/operator-control.json" > "$PUNCTURE_RUN/operator-control.log" 2>&1
python3 "$PUNCTURE_ROOT/workflow-axisym/validation/run_trumpet_pilot.py" --library "$PUNCTURE_ROOT/build-radial512-cuda/libHiSpID.so" --case gamma10 --axisymmetric --n 512 --npolar 320 --nphi 8 --krylov-restart 64 --memory-mib 65536 --initial "$PUNCTURE_ROOT/axisymmetric/gamma10_256x320/diagnostic.checkpoint" --initial-source-library "$PUNCTURE_ROOT/build-axisym-cuda/libHiSpID.so" --output "$PUNCTURE_RUN/gamma10_512x320" > "$PUNCTURE_RUN/solve.log" 2>&1
