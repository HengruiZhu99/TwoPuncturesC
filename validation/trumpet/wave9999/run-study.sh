#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
cd "$base/validation-wave9999-v1"
export OMP_NUM_THREADS=16 OMP_PROC_BIND=spread OMP_PLACES=cores OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
python validation/wave22_study.py --library "$base/build-wave22-v1/libHiSpID.so" --output "$base/wave9999-study-v1" --levels 96 128 160 --focus-separation 32 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers > "$base/wave9999-study-v1.log" 2>&1
