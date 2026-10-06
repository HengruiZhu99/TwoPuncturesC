#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/analytic-matrices"
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source-analytic/python"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_ROOT/source-analytic" -B "$PUNCTURE_ROOT/build-analytic-cuda" -DCMAKE_BUILD_TYPE=Release -DHISPID_ANALYTIC_MATRICES=ON -DHISPID_ROW_POWER=3 -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper > "$PUNCTURE_RUN/configure.log" 2>&1
/usr/bin/time -v cmake --build "$PUNCTURE_ROOT/build-analytic-cuda" --target HiSpID test_hispid_axis -j4 > "$PUNCTURE_RUN/build.log" 2>&1
"$PUNCTURE_ROOT/build-analytic-cuda/test_hispid_axis" > "$PUNCTURE_RUN/axis.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-analytic-cuda/libHiSpID.so" "$PUNCTURE_ROOT/build-analytic-cuda/libTwoPunctures.so" "$PUNCTURE_ROOT/source-analytic/src/HiSpID_spectral.hpp" "$PUNCTURE_ROOT/source-analytic/src/HiSpID_setup_kokkos.cpp" "$PUNCTURE_ROOT/source-analytic/CMakeLists.txt" > "$PUNCTURE_RUN/images-sources.sha256"
python3 "$PUNCTURE_ROOT/workflow-polar/validation/check_trumpet_production.py" --library "$PUNCTURE_ROOT/build-analytic-cuda/libHiSpID.so" --geometry host execution --skip-seeds --output "$PUNCTURE_RUN/control.json" > "$PUNCTURE_RUN/control.log" 2>&1
python3 -c 'import json,sys;sys.exit(not json.load(open(sys.argv[1]))["passed"])' "$PUNCTURE_RUN/control.json"
python3 "$PUNCTURE_ROOT/workflow-polar/validation/run_trumpet_tolerance.py" --library "$PUNCTURE_ROOT/build-analytic-cuda/libHiSpID.so" --n 256 --npolar 512 --nphi 16 --krylov-restart 64 --memory-mib 65536 --initial "$PUNCTURE_ROOT/convergence/moderate256/diagnostic.checkpoint" --initial-source-library "$PUNCTURE_ROOT/build-polar-cuda/libHiSpID.so" --output "$PUNCTURE_RUN/moderate256" > "$PUNCTURE_RUN/solve.log" 2>&1
