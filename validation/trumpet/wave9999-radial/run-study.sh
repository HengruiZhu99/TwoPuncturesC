#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
cd "$base/validation-wave-radial-map-v2"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
HISPID_TRACE_NEWTON=1 python validation/wave22_study.py --library "$base/build-wave-radial-map-v1/libHiSpID.so" --output "$base/wave9999-radial-map-study-v6" --levels 384 480 576 --polar-points 128 --memory-limit-mib 65536 --focus-separation 16 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers --linear-rtol .02 --max-krylov 4000 --max-newton 40 --krylov-restart 200 --krylov gmres --polish-linear-rtol .5 --initial-checkpoint "$base/wave9999-radial-map-study-v2/n384.hispid" > "$base/wave9999-radial-map-study-v6.log" 2>&1
