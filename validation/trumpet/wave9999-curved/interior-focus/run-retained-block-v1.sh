#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
flags=(-DCMAKE_BUILD_TYPE=Release -DHISPID_AXIS_TAU=ON -DHISPID_FOCUS_TAU=ON -DHISPID_FOCUS_ORDER=1 -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_STREAM_POLAR_PRECONDITIONER=ON -DHISPID_ROW_POWER=3 -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_BATCHED_MODAL_TRANSFER=ON -DHISPID_RADIAL_STRETCH=16 -DHISPID_ANGULAR_STRETCH=.3 -DHISPID_CURVED_MODAL_PRECONDITIONER=ON -DHISPID_CONSTRAINT_BLOCK_PRECONDITIONER=ON)
cmake -S "$base/source-wave-curved-v16" -B "$base/build-wave-curved-v16" "${flags[@]}" -DPUNCTURES_KOKKOS=ON -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper -DKokkos_DIR="$base/build-cuda/cmake_packages/Kokkos"
cmake --build "$base/build-wave-curved-v16" --target HiSpID test_kokkos test_hispid_constraint_block -j16
/opt/cray/pe/gcc-native/14/bin/gcc -std=c99 -O2 -I"$base/source-wave-curved-v16/include" "$base/source-wave-curved-v16/src/PunctureKrylov.c" "$base/source-wave-curved-v16/tests/test_krylov.c" -lm -o "$base/build-wave-curved-v16/test_retained_krylov"
"$base/build-wave-curved-v16/test_retained_krylov"
"$base/build-wave-curved-v16/test_kokkos"
HISPID_TRACE_NEWTON=1 "$base/build-wave-curved-v16/test_hispid_constraint_block" 1
sha256sum "$base/build-wave-curved-v16/libHiSpID.so" "$base/build-wave-curved-v16/libTwoPunctures.so"
cd "$base/source-wave-curved-v16"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
set +e
HISPID_TRACE_NEWTON=1 python validation/wave22_study.py --library "$base/build-wave-curved-v16/libHiSpID.so" --output "$base/wave9999-retained-block-seed-v1" --levels 384 --polar-levels 128 --memory-limit-mib 98304 --focus-separation 1 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers --linear-rtol .1 --max-krylov 40 --max-newton 2 --krylov-restart 12 --krylov gmres --initial-checkpoint "$base/wave9999-toroidal-seed-v1/n384.hispid"
set -e
python - "$base/wave9999-retained-block-seed-v1/result.json" <<'PYSEED'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);d=json.loads(p.read_text());assert len(d['levels'])==1
r=d['levels'][0];assert (p.parent/'n384.hispid').is_file()
assert max(r['diagnostics']['scaled_linf'])<r['initial_scaled_linf'],r['diagnostics']
assert r['diagnostics']['converged'] or r['diagnostics']['error']=='Newton iteration limit',r['diagnostics']
print('retained fixed-basis diagnostic progress',r['diagnostics'],flush=True)
PYSEED
HISPID_TRACE_NEWTON=1 python validation/wave22_study.py --library "$base/build-wave-curved-v16/libHiSpID.so" --output "$base/wave9999-retained-block-study-v1" --levels 384 480 576 --polar-levels 128 160 192 --memory-limit-mib 98304 --focus-separation 1 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers --linear-rtol .1 --max-krylov 300 --max-newton 32 --krylov-restart 12 --krylov gmres --initial-checkpoint "$base/wave9999-retained-block-seed-v1/n384.hispid"
