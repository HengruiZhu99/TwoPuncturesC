#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s  %s\n' 1eec375773a0fd8f287320b08972364e075441c35fb4b1c1623d48ceeff1a4b0 "$base/build-wave-curved-v17/libHiSpID.so" | sha256sum -c -
cd "$base/study-wave-curved-v17"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
HISPID_TRACE_NEWTON=1 python validation/wave22_study.py --library "$base/build-wave-curved-v17/libHiSpID.so" --output "$base/wave9999-inner-augmentation-study-v1" --levels 384 480 576 --polar-levels 128 160 192 --memory-limit-mib 98304 --focus-separation 1 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers --linear-rtol .1 --max-krylov 80 --max-newton 32 --krylov-restart 24 --krylov gmres --checkpoint-each-newton --prolongate --initial-checkpoint "$base/wave9999-inner-augmentation-pilot-v1/n384.hispid"
bash "$base/run-inner-augmentation-post-v1.sh" "$base/wave9999-inner-augmentation-study-v1" "$base/build-wave-curved-v17/libHiSpID.so" wave9999-inner-augmentation-final-v1
