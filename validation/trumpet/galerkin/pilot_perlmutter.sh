#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE="$PUNCTURE_ROOT/source-galerkin-20261018"
PUNCTURE_BUILD="$PUNCTURE_ROOT/build-galerkin-20261018"
PUNCTURE_RUN="$PUNCTURE_ROOT/galerkin-pilot-v2"
mkdir "$PUNCTURE_RUN"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_SOURCE" -B "$PUNCTURE_BUILD" \
 -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_REGULARITY_CAP=4 \
 -DHISPID_AXIS_TAU=ON -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_COMPENSATED_MODAL_PROJECTION=ON \
 -DHISPID_OVERINTEGRATED_GALERKIN=ON -DHISPID_BATCHED_MODAL_TRANSFER=ON \
 -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" \
 -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc \
 -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_BUILD" --target HiSpID test_hispid_galerkin_transfer -j4 > "$PUNCTURE_RUN/build.log" 2>&1
sha256sum "$PUNCTURE_BUILD/libHiSpID.so" "$PUNCTURE_BUILD/libTwoPunctures.so" > "$PUNCTURE_RUN/images.sha256"
"$PUNCTURE_BUILD/test_hispid_galerkin_transfer" > "$PUNCTURE_RUN/transfer-control.txt"
python3 "$PUNCTURE_SOURCE/validation/trumpet/c4-experiment/scaled-inverse/control.py" \
 --library "$PUNCTURE_BUILD/libHiSpID.so" --execution kokkos --solve \
 --output "$PUNCTURE_RUN/control.json" > "$PUNCTURE_RUN/control.log" 2>&1
python3 - "$PUNCTURE_RUN/control.json" "$PUNCTURE_SOURCE/validation/trumpet/galerkin/control/cpu-unweighted.json" <<'CHECK'
import json,sys,numpy as np
from pathlib import Path
r=json.load(open(sys.argv[1]));c=json.load(open(sys.argv[2]));assert r['diagnostics']['converged'] and c['diagnostics']['converged']
e=max(float(np.max(abs(np.array(r['field_witness'][k])-np.array(c['field_witness'][k]))/(1+abs(np.array(c['field_witness'][k]))))) for k in c['field_witness'])
Path(sys.argv[1]).with_name('backend-comparison.json').write_text(json.dumps(dict(field_scaled_difference=e,limit=1e-10))+'\n');assert e<1e-10
CHECK
export HISPID_TRACE_NEWTON=1
python3 "$PUNCTURE_SOURCE/validation/run_trumpet_pilot.py" \
 --library "$PUNCTURE_BUILD/libHiSpID.so" --case moderate \
 --n 192 --npolar 384 --nphi 16 --krylov-restart 32 \
 --memory-mib 65536 --linear-rtol .1 \
 --initial "$PUNCTURE_ROOT/axis-tau-polar-full-v1/moderate192/diagnostic.checkpoint" \
 --initial-source-library "$PUNCTURE_ROOT/build-axis-tau-polar/libHiSpID.so" \
 --output "$PUNCTURE_RUN/moderate192" > "$PUNCTURE_RUN/solve192.log" 2>&1
