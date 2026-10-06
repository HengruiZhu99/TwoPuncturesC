#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export PYTHONPATH="$PUNCTURE_ROOT/source-batched-transfer/python:$PUNCTURE_ROOT/source-batched-transfer/validation:$PUNCTURE_ROOT/source-batched-transfer/examples"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python3 "$PUNCTURE_ROOT/source-batched-transfer/validation/check_trumpet_axisymmetric.py" --library "$PUNCTURE_ROOT/build-batched-transfer-cuda/libHiSpID.so" --execution kokkos --output "$PUNCTURE_ROOT/batched-transfer-cuda/control-v2.json" > "$PUNCTURE_ROOT/batched-transfer-cuda/control-v2.log" 2>&1
# Import the queued C4 pilot too, before spending an allocation on it.
PYTHONPATH="$PUNCTURE_ROOT/source-c4-experiment/python:$PUNCTURE_ROOT/source-c4-experiment/validation:$PUNCTURE_ROOT/source-c4-experiment/examples" python3 "$PUNCTURE_ROOT/source-c4-experiment/validation/run_trumpet_pilot.py" --help > "$PUNCTURE_ROOT/c4-experiment/pilot-import.txt"
