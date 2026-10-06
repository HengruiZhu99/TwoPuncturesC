#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/axis-tau-spectral-radial-v1"
mkdir "$PUNCTURE_RUN"
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_ROOT/source-axis-tau-spectral-radial" -B "$PUNCTURE_ROOT/build-axis-tau-spectral-radial" \
  -DCMAKE_BUILD_TYPE=Release -DHISPID_ROW_POWER=3 -DHISPID_REGULARITY_CAP=4 \
  -DHISPID_RADIAL_STRETCH=.2 -DHISPID_ANGULAR_STRETCH=2. -DHISPID_NEWTON_BACKTRACKS=10 -DHISPID_STABLE_SCALAR_SOURCE=OFF -DHISPID_MONOTONE_PRECONDITIONER=OFF -DHISPID_AXIS_TAU=ON -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_SPECTRAL_RADIAL_PRECONDITIONER=ON \
  -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" \
  -DCMAKE_C_COMPILER=/opt/cray/pe/gcc-native/14/bin/gcc \
  -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-axis-tau-spectral-radial" --target HiSpID test_hispid_radial_preconditioner -j4 > "$PUNCTURE_RUN/build.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-axis-tau-spectral-radial/libHiSpID.so" "$PUNCTURE_ROOT/source-axis-tau-spectral-radial/src/HiSpID_solver.cpp" > "$PUNCTURE_RUN/images-sources.sha256"
"$PUNCTURE_ROOT/build-axis-tau-spectral-radial/test_hispid_radial_preconditioner" > "$PUNCTURE_RUN/polynomial-control.txt"
python3 "$PUNCTURE_ROOT/source-axis-tau-spectral-radial/validation/trumpet/c4-experiment/scaled-inverse/control.py" \
 --library "$PUNCTURE_ROOT/build-axis-tau-spectral-radial/libHiSpID.so" --execution kokkos \
 --output "$PUNCTURE_RUN/control.json" --solve --krylov gmres > "$PUNCTURE_RUN/control.log" 2>&1
python3 - "$PUNCTURE_RUN/control.json" "$PUNCTURE_ROOT/source-axis-tau-spectral-radial/validation/trumpet/c4-experiment/axis-tau/polar-exact/spectral-radial/cached-control.json" <<'CHECK'
import json,sys,numpy as np
from pathlib import Path
r=json.load(open(sys.argv[1]));assert r['diagnostics']['converged'] and r['resolved_options']['krylov']=='gmres'
c=json.load(open(sys.argv[2]));assert c['diagnostics']['converged']
error=max(float(np.max(abs(np.array(r['field_witness'][k])-np.array(c['field_witness'][k]))/(1+abs(np.array(c['field_witness'][k]))))) for k in c['field_witness'])
Path(sys.argv[1]).with_name('backend-comparison.json').write_text(json.dumps(dict(field_scaled_difference=error,limit=1e-10))+'\n')
assert error<1e-10
CHECK
python3 "$PUNCTURE_ROOT/source-axis-tau-spectral-radial/validation/diagnose_trumpet_newton.py" \
 --library "$PUNCTURE_ROOT/build-axis-tau-spectral-radial/libHiSpID.so" \
 --source-library "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" \
 --checkpoint "$PUNCTURE_ROOT/axis-tau-compensated-v3/moderate240/diagnostic.checkpoint" \
 --krylov gmres --linear-rtol .1 --max-krylov 2400 --krylov-restart 64 \
 --output "$PUNCTURE_RUN/probe.json" > "$PUNCTURE_RUN/probe.log" 2>&1
