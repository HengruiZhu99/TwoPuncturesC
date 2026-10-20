#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
expected=1eec375773a0fd8f287320b08972364e075441c35fb4b1c1623d48ceeff1a4b0
actual=$(sha256sum "$base/build-wave-curved-v17/libHiSpID.so" | cut -d' ' -f1)
test "$actual" = "$expected"
cd "$base/source-wave-curved-v17"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
export PUNCTURE_KRYLOV_ARNOLDI="$base/frozen-arnoldi-v1.jsonl"
set +e
HISPID_TRACE_NEWTON=1 HISPID_PROBE_LINEAR_ACTION=1 python validation/wave22_study.py --library "$base/build-wave-curved-v17/libHiSpID.so" --output "$base/wave9999-frozen-arnoldi-v1" --levels 384 --polar-levels 128 --memory-limit-mib 98304 --focus-separation 1 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers --linear-rtol .1 --max-krylov 1 --max-newton 1 --krylov-restart 24 --krylov gmres --initial-checkpoint "$base/wave9999-inner-augmentation-pilot-v1/n384.hispid"
status=$?
printf 'bounded_physical_driver_exit=%s\n' "$status"
exit "$status"
