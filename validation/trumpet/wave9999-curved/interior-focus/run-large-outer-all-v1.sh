#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s  %s\n' fa97396c525bc8663d2b65d0fa80320a2aa7f68c5dcf8627bcd28fbe32935e54 "$base/build-wave-curved-v16/libHiSpID.so" | sha256sum -c -
cd "$base/source-wave-curved-v16"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
HISPID_TRACE_NEWTON=1 python validation/wave22_study.py --library "$base/build-wave-curved-v16/libHiSpID.so" --output "$base/wave9999-retained-large-outer-study-v1" --levels 384 480 576 --polar-levels 128 160 192 --memory-limit-mib 98304 --focus-separation 1 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers --linear-rtol .1 --max-krylov 300 --max-newton 32 --krylov-restart 24 --krylov lgmres --initial-checkpoint "$base/wave9999-retained-block-seed-v1/n384.hispid"
bash "$base/run-retained-post-v1.sh" "$base/wave9999-retained-large-outer-study-v1" "$base/build-wave-curved-v16/libHiSpID.so" wave9999-retained-large-outer-final
