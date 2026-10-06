#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/focused-map"
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source/python"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_ROOT/source-focused" -B "$PUNCTURE_ROOT/build-focused" -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_RADIAL_STRETCH=.05 -DHISPID_ANGULAR_STRETCH=3. -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc > "$PUNCTURE_RUN/configure.log" 2>&1
/usr/bin/time -v cmake --build "$PUNCTURE_ROOT/build-focused" --target HiSpID -j4 > "$PUNCTURE_RUN/build.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-focused/libHiSpID.so" "$PUNCTURE_ROOT/build-focused/libTwoPunctures.so" > "$PUNCTURE_RUN/images.sha256"
python3 "$PUNCTURE_ROOT/workflow-focused/validation/check_trumpet_production.py" --library "$PUNCTURE_ROOT/build-focused/libHiSpID.so" --geometry host --skip-seeds --output "$PUNCTURE_RUN/control.json" > "$PUNCTURE_RUN/control.log" 2>&1
python3 -c 'import json,sys;sys.exit(not json.load(open(sys.argv[1]))["passed"])' "$PUNCTURE_RUN/control.json"
python3 "$PUNCTURE_ROOT/workflow-focused/validation/run_trumpet_pilot.py" --library "$PUNCTURE_ROOT/build-focused/libHiSpID.so" --n 160 --npolar 256 --nphi 16 --memory-mib 32768 --initial "$PUNCTURE_ROOT/refinement-physics/moderate256/diagnostic.checkpoint" --initial-source-library "$PUNCTURE_ROOT/build-cuda/libHiSpID.so" --output "$PUNCTURE_RUN/moderate160x256" > "$PUNCTURE_RUN/solve.log" 2>&1
