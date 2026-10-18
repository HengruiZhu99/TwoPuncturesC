#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE="$PUNCTURE_ROOT/source-arnoldi-20261018"
PUNCTURE_BUILD="$PUNCTURE_ROOT/build-arnoldi-20261018"
PUNCTURE_RUN="$PUNCTURE_ROOT/arnoldi-retained240-v1"
mkdir "$PUNCTURE_RUN"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_SOURCE" -B "$PUNCTURE_BUILD" \
 -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_REGULARITY_CAP=4 \
 -DHISPID_AXIS_TAU=ON -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_COMPENSATED_MODAL_PROJECTION=ON \
 -DHISPID_OVERINTEGRATED_GALERKIN=OFF -DHISPID_DAMPED_MODAL_RELAXATION=OFF \
 -DHISPID_BATCHED_MODAL_TRANSFER=ON -DPUNCTURES_KOKKOS=ON \
 -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" \
 -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc \
 -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_BUILD" --target HiSpID -j4 > "$PUNCTURE_RUN/build.log" 2>&1
sha256sum "$PUNCTURE_BUILD/libHiSpID.so" "$PUNCTURE_BUILD/libTwoPunctures.so" > "$PUNCTURE_RUN/images.sha256"
export PUNCTURE_KRYLOV_ARNOLDI="$PUNCTURE_RUN/arnoldi.jsonl"
python3 "$PUNCTURE_SOURCE/validation/diagnose_trumpet_newton.py" \
 --library "$PUNCTURE_BUILD/libHiSpID.so" \
 --source-library "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" \
 --checkpoint "$PUNCTURE_ROOT/axis-tau-compensated-v3/moderate240/diagnostic.checkpoint" \
 --krylov gmres --linear-rtol .1 --max-krylov 128 --krylov-restart 128 --save-step \
 --output "$PUNCTURE_RUN/probe.json" > "$PUNCTURE_RUN/probe.log" 2>&1
