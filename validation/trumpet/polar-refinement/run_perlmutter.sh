#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/polar-refinement"
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source-polar/python"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
for PUNCTURE_STAGE in sampler cuda; do
  PUNCTURE_FLAGS=(-DPUNCTURES_KOKKOS=OFF -DCMAKE_CXX_COMPILER="$NVCC_WRAPPER_DEFAULT_COMPILER")
  if [[ "$PUNCTURE_STAGE" == cuda ]]; then
    PUNCTURE_FLAGS=(-DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper)
  fi
  cmake -S "$PUNCTURE_ROOT/source-polar" -B "$PUNCTURE_ROOT/build-polar-$PUNCTURE_STAGE" -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_RADIAL_STRETCH=.2 -DHISPID_ANGULAR_STRETCH=2. -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc "${PUNCTURE_FLAGS[@]}" > "$PUNCTURE_RUN/configure-$PUNCTURE_STAGE.log" 2>&1
  /usr/bin/time -v cmake --build "$PUNCTURE_ROOT/build-polar-$PUNCTURE_STAGE" --target HiSpID -j4 > "$PUNCTURE_RUN/build-$PUNCTURE_STAGE.log" 2>&1
  sha256sum "$PUNCTURE_ROOT/build-polar-$PUNCTURE_STAGE/libHiSpID.so" "$PUNCTURE_ROOT/build-polar-$PUNCTURE_STAGE/libTwoPunctures.so" > "$PUNCTURE_RUN/images-$PUNCTURE_STAGE.sha256"
done
python3 "$PUNCTURE_ROOT/athenak-polar/tst/test_suite/z4c/check_trumpet_checkpoint.py" --native-root "$PUNCTURE_ROOT/source-polar" --library "$PUNCTURE_ROOT/build-polar-sampler/libHiSpID.so" --polar-extent 384 --output "$PUNCTURE_RUN/reader384" > "$PUNCTURE_RUN/reader.log" 2>&1
python3 "$PUNCTURE_ROOT/workflow-polar/validation/run_trumpet_pilot.py" --library "$PUNCTURE_ROOT/build-polar-cuda/libHiSpID.so" --n 256 --npolar 384 --nphi 16 --memory-mib 65536 --initial "$PUNCTURE_ROOT/refinement-physics/moderate256/diagnostic.checkpoint" --initial-source-library "$PUNCTURE_ROOT/build-cuda/libHiSpID.so" --output "$PUNCTURE_RUN/moderate256x384" > "$PUNCTURE_RUN/solve.log" 2>&1
