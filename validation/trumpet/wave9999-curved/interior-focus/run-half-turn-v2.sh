#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s  %s\n' b9517e9895f5f0a495e1656d73f2572a0806a6b4742705fbc5b983a715672bb6 "$base/build-wave-curved-v26/libHiSpID.so" | sha256sum -c -
cd "$base/study-wave-curved-v17"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
# Cold start: original saved nonsymmetric diagnostic is deliberately not loaded.
set +e
HISPID_TRACE_NEWTON=1 python validation/wave22_study.py --library "$base/build-wave-curved-v26/libHiSpID.so" --output "$base/wave9999-half-turn-pilot-v1" --levels 384 --polar-levels 128 --memory-limit-mib 98304 --focus-separation 1 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers --linear-rtol .1 --max-krylov 6 --max-newton 8 --krylov-restart 24 --krylov gmres --checkpoint-each-newton
status=$?
printf 'bounded_physical_driver_exit=%s\n' "$status"
