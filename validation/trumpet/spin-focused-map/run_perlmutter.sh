#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/spin-focused-map"
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source-polar/python"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_ROOT/source-polar" -B "$PUNCTURE_ROOT/build-spin-focused" -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_RADIAL_STRETCH=.03 -DHISPID_ANGULAR_STRETCH=3.5 -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper > "$PUNCTURE_RUN/configure.log" 2>&1
/usr/bin/time -v cmake --build "$PUNCTURE_ROOT/build-spin-focused" --target HiSpID -j4 > "$PUNCTURE_RUN/build.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-spin-focused/libHiSpID.so" "$PUNCTURE_ROOT/build-spin-focused/libTwoPunctures.so" "$PUNCTURE_ROOT/source-polar/CMakeLists.txt" > "$PUNCTURE_RUN/images-sources.sha256"
python3 "$PUNCTURE_ROOT/workflow-polar/validation/check_trumpet_production.py" --library "$PUNCTURE_ROOT/build-spin-focused/libHiSpID.so" --geometry host --skip-seeds --output "$PUNCTURE_RUN/control.json" > "$PUNCTURE_RUN/control.log" 2>&1
python3 -c 'import json,sys;sys.exit(not json.load(open(sys.argv[1]))["passed"])' "$PUNCTURE_RUN/control.json"
python3 "$PUNCTURE_ROOT/workflow-polar/validation/run_trumpet_tolerance.py" --library "$PUNCTURE_ROOT/build-spin-focused/libHiSpID.so" --case spin99 --n 160 --npolar 320 --nphi 32 --krylov-restart 64 --memory-mib 65536 --initial "$PUNCTURE_ROOT/extreme-pilots/spin99_160/diagnostic.checkpoint" --initial-source-library "$PUNCTURE_ROOT/build-polar-cuda/libHiSpID.so" --output "$PUNCTURE_RUN/spin99_160" > "$PUNCTURE_RUN/solve.log" 2>&1
