#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s  %s\n' ddc17732c4b1f88f86acc8bc01ac6147d62cb5a2f6b2ca7450f97fef5443f38c "$base/build-wave-curved-v27/libHiSpID.so" | sha256sum -c -
cd "$base/study-wave-curved-v17"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
HISPID_TRACE_NEWTON=1 python validation/wave22_study.py --library "$base/build-wave-curved-v27/libHiSpID.so" --output "$base/wave9999-half-turn-cold-study-v1" --levels 384 480 576 --polar-levels 128 160 192 --memory-limit-mib 98304 --focus-separation 1 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers --linear-rtol .1 --max-krylov 6 --max-newton 12 --krylov-restart 24 --krylov gmres --checkpoint-each-newton --initial-checkpoint "$base/wave9999-half-turn-pilot-v1/n384.hispid"
bash "$base/run-inner-augmentation-post-v1.sh" "$base/wave9999-half-turn-cold-study-v1" "$base/build-wave-curved-v27/libHiSpID.so" wave9999-half-turn-cold-final-v1
