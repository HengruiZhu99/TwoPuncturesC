#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
cd "$base/source-wave-radial-map-v1"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python validation/wave22_radial_profiles.py --library "$base/build-wave-radial-map-cpu-v1/libHiSpID.so" --self-test --output "$base/wave9999-radial-map-oracle-v1.json" > "$base/wave9999-radial-map-oracle-v1.log" 2>&1
python validation/wave22_study.py --library "$base/build-wave-radial-map-v1/libHiSpID.so" --output "$base/wave9999-radial-map-study-v1" --levels 384 480 576 --polar-points 128 --memory-limit-mib 65536 --focus-separation 16 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers --linear-rtol .02 --max-krylov 300 --polish-linear-rtol .5 > "$base/wave9999-radial-map-study-v1.log" 2>&1
