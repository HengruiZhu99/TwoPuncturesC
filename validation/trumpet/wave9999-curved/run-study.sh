#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
"$base/build-wave-curved-v3/test_by_kokkos"
sha256sum "$base/build-wave-curved-v3/libHiSpID.so"
cd "$base/source-wave-curved-v3"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
HISPID_TRACE_NEWTON=1 python validation/wave22_study.py --library "$base/build-wave-curved-v3/libHiSpID.so" --output "$base/wave9999-curved-study-v3" --levels 384 480 576 --polar-points 128 --memory-limit-mib 65536 --focus-separation 16 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers --linear-rtol .1 --max-krylov 1200 --max-newton 32 --krylov-restart 100 --krylov gmres --initial-checkpoint "$base/wave9999-radial-map-study-v2/n384.hispid"
