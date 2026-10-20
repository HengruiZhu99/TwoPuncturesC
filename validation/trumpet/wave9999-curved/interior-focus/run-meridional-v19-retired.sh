#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
flags=(-DCMAKE_BUILD_TYPE=Release -DHISPID_AXIS_TAU=ON -DHISPID_FOCUS_TAU=ON -DHISPID_FOCUS_ORDER=1 -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_STREAM_POLAR_PRECONDITIONER=ON -DHISPID_ROW_POWER=3 -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_BATCHED_MODAL_TRANSFER=ON -DHISPID_RADIAL_STRETCH=16 -DHISPID_ANGULAR_STRETCH=.3 -DHISPID_CURVED_MODAL_PRECONDITIONER=ON -DHISPID_CONSTRAINT_BLOCK_PRECONDITIONER=ON -DHISPID_MERIDIONAL_PRECONDITIONER=ON)
cmake -S "$base/source-wave-curved-v19" -B "$base/build-wave-curved-v19" "${flags[@]}" -DPUNCTURES_KOKKOS=ON -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper -DKokkos_DIR="$base/build-cuda/cmake_packages/Kokkos"
cmake --build "$base/build-wave-curved-v19" --target test_hispid_meridional test_hispid_constraint_block -j16
"$base/build-wave-curved-v19/test_hispid_meridional"
"$base/build-wave-curved-v19/test_hispid_constraint_block" 1
sha256sum "$base/build-wave-curved-v19/libHiSpID.so"
cd "$base/source-wave-curved-v19"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
# Diagnostic includes actual frozen physical RHS; the existing smooth control
# is only an algebraic gate. Failed native acceptance remains explicit.
set +e
HISPID_TRACE_NEWTON=1 HISPID_PROBE_LINEAR_ACTION=1 python validation/wave22_study.py --library "$base/build-wave-curved-v19/libHiSpID.so" --output "$base/wave9999-meridional-pilot-v1" --levels 384 --polar-levels 128 --memory-limit-mib 98304 --focus-separation 1 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers --linear-rtol .1 --max-krylov 6 --max-newton 1 --krylov-restart 24 --krylov gmres --initial-checkpoint "$base/wave9999-inner-augmentation-pilot-v1/n384.hispid"
status=$?
printf 'bounded_physical_driver_exit=%s\n' "$status"
# A bounded Newton-limit result is usable only as an initial guess after
# an actually accepted update; no convergence status is transferred.
python - "$base/wave9999-meridional-pilot-v1/result.json" <<'PYGATE'
import json,sys
r=json.load(open(sys.argv[1]));d=r['levels'][0]['diagnostics']
accepted=bool(d['converged']) or (d['error']=='Newton iteration limit' and d['newton_iterations']==1 and max(d['scaled_linf'])<r['levels'][0]['initial_scaled_linf'])
if not accepted:raise SystemExit('bounded physical pilot did not accept a Newton update; full study not started')
PYGATE
gate=$?
set -e
if [ "$gate" -ne 0 ]; then exit "$gate"; fi
bash "$base/run-meridional-study-v1.sh"

