#!/usr/bin/env bash
# Run through srun inside a one-GPU shared_interactive allocation. Every
# configure/build/test is sequential; numerical tests use one host thread.
set -euo pipefail
export PUNCTURE_RUN_ROOT=${1:?isolated scratch run directory required}
cd "$PUNCTURE_RUN_ROOT/source"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores
export PYTHONPATH=python:validation:examples
if [[ ! -d "$PUNCTURE_RUN_ROOT/frozen-reference" ]]; then
  python3 validation/prepare_kokkos_reference.py --archive "$PUNCTURE_RUN_ROOT/frozen-25ca064.tar" --destination "$PUNCTURE_RUN_ROOT/frozen-reference" > "$PUNCTURE_RUN_ROOT/reference-provenance.log"
fi
cp validation/by_solver_timer.c "$PUNCTURE_RUN_ROOT/frozen-reference/validation/by_solver_timer.c"
python3 - <<'PY'
import hashlib,json,os
from pathlib import Path
root=Path(os.environ['PUNCTURE_RUN_ROOT'])/'frozen-reference';path=root/'reference-build-provenance.json';record=json.loads(path.read_text())
record['timing_wrapper_sha256']=hashlib.sha256((root/'validation/by_solver_timer.c').read_bytes()).hexdigest()
record['timing_note']='Current wrappers add native residual timing and preserve reference arithmetic.'
path.write_text(json.dumps(record,indent=2)+'\n')
PY
cmake -S "$PUNCTURE_RUN_ROOT/frozen-reference" -B "$PUNCTURE_RUN_ROOT/build-reference" -DCMAKE_BUILD_TYPE=Release -DPUNCTURES_KOKKOS=OFF -DHISPID_ROW_POWER=3 -DPUNCTURES_BENCHMARK=ON -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ > "$PUNCTURE_RUN_ROOT/configure-reference.log" 2>&1
cmake --build "$PUNCTURE_RUN_ROOT/build-reference" -j1 > "$PUNCTURE_RUN_ROOT/build-reference.log" 2>&1
for PUNCTURE_SPACE in serial openmp; do
  PUNCTURE_OMP=OFF
  if [[ "$PUNCTURE_SPACE" == openmp ]]; then PUNCTURE_OMP=ON; fi
  cmake -S . -B "$PUNCTURE_RUN_ROOT/build-$PUNCTURE_SPACE" -DCMAKE_BUILD_TYPE=Release -DPUNCTURES_KOKKOS=ON -DPUNCTURES_KOKKOS_SOURCE="$PUNCTURE_RUN_ROOT/kokkos" -DHISPID_ROW_POWER=3 -DPUNCTURES_BENCHMARK=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP="$PUNCTURE_OMP" -DKokkos_ENABLE_CUDA=OFF -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ > "$PUNCTURE_RUN_ROOT/configure-$PUNCTURE_SPACE.log" 2>&1
  cmake --build "$PUNCTURE_RUN_ROOT/build-$PUNCTURE_SPACE" -j1 > "$PUNCTURE_RUN_ROOT/build-$PUNCTURE_SPACE.log" 2>&1
  ctest --test-dir "$PUNCTURE_RUN_ROOT/build-$PUNCTURE_SPACE" --output-on-failure -j1 > "$PUNCTURE_RUN_ROOT/test-$PUNCTURE_SPACE.log" 2>&1
 done
cmake --build "$PUNCTURE_RUN_ROOT/build-cuda" -j1 > "$PUNCTURE_RUN_ROOT/build-cuda.log" 2>&1
ctest --test-dir "$PUNCTURE_RUN_ROOT/build-cuda" --output-on-failure -j1 > "$PUNCTURE_RUN_ROOT/test-cuda.log" 2>&1
python3 validation/kokkos_build_manifest.py --root "$PUNCTURE_RUN_ROOT" --output "$PUNCTURE_RUN_ROOT/build-manifest.json"
