#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
flags=(-DCMAKE_BUILD_TYPE=Release -DHISPID_AXIS_TAU=ON -DHISPID_FOCUS_TAU=ON -DHISPID_FOCUS_ORDER=1 -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_STREAM_POLAR_PRECONDITIONER=ON -DHISPID_ROW_POWER=3 -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_BATCHED_MODAL_TRANSFER=ON -DHISPID_RADIAL_STRETCH=16 -DHISPID_ANGULAR_STRETCH=.3 -DHISPID_CURVED_MODAL_PRECONDITIONER=ON -DHISPID_CONSTRAINT_BLOCK_PRECONDITIONER=ON -DHISPID_SPECTRAL_POLAR_PRECONDITIONER=ON)
cmake -S "$base/source-wave-curved-v21" -B "$base/build-wave-curved-v21" "${flags[@]}" -DPUNCTURES_KOKKOS=ON -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper -DKokkos_DIR="$base/build-cuda/cmake_packages/Kokkos"
cmake --build "$base/build-wave-curved-v21" --target HiSpID test_hispid_spectral_polar test_hispid_constraint_block -j16
sha256sum "$base/build-wave-curved-v21/libHiSpID.so"
"$base/build-wave-curved-v21/test_hispid_spectral_polar" 4
"$base/build-wave-curved-v21/test_hispid_spectral_polar" 16
"$base/build-wave-curved-v21/test_hispid_constraint_block" 1
cd "$base/source-wave-curved-v21"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
args=(--library "$base/build-wave-curved-v21/libHiSpID.so" --levels 384 --polar-levels 128 --memory-limit-mib 98304 --focus-separation 1 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers --linear-rtol .1 --max-newton 1 --krylov-restart 24 --krylov gmres --initial-checkpoint "$base/wave9999-inner-augmentation-pilot-v1/n384.hispid")
set +e
HISPID_TRACE_NEWTON=1 HISPID_PROBE_MODAL_RESPONSE=1 python validation/wave22_study.py "${args[@]}" --max-krylov 1 --output "$base/wave9999-spectral-polar-audit-v1" > "$base/spectral-polar-audit-v1.log" 2>&1
status=$?
set -e
cat "$base/spectral-polar-audit-v1.log"
printf 'diagnostic_driver_exit=%s\n' "$status"
python - "$base" <<'CHECK'
import json,re,sys
from pathlib import Path
base=Path(sys.argv[1]);d=json.loads((base/'wave9999-spectral-polar-audit-v1/result.json').read_text())['levels'][0]['diagnostics']
assert not d['converged'] and d['krylov_iterations']==0 and d['error']=='modal response diagnostic only; no Krylov solve attempted'
log=(base/'spectral-polar-audit-v1.log').read_text();m=re.search(r'sum_action_l2=([\deE.+-]+).*superposition_relative=([\deE.+-]+)',log);assert m
ratio=float(m[1])/143653246.78657773
passed=ratio<.01 and float(m[2])<1e-10
receipt=dict(passed=passed,combined_action_ratio=ratio,scope='Approximate-inverse diagnostic gate only; not physical acceptance',next='one bounded Newton attempt' if passed else 'stop: amplification not substantially improved')
(base/'spectral-polar-audit-gate-v1.json').write_text(json.dumps(receipt,indent=2)+'\n');print(receipt,flush=True)
CHECK
if ! python - "$base/spectral-polar-audit-gate-v1.json" <<'CHECK'
import json,sys
sys.exit(0 if json.load(open(sys.argv[1]))['passed'] else 1)
CHECK
then
 echo 'Candidate stopped at diagnostic gate; no physical solve attempted.'
 exit 0
fi
set +e
HISPID_TRACE_NEWTON=1 python validation/wave22_study.py "${args[@]}" --max-krylov 6 --output "$base/wave9999-spectral-polar-pilot-v1"
status=$?
set -e
printf 'bounded_physical_driver_exit=%s\n' "$status"
echo 'Bounded pilot finished. Three-grid continuation requires review of accepted update and finest-grid memory plan.'
