#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/affine-horizon/exact"
mkdir -p "$PUNCTURE_RUN"
export PYTHONPATH="$PUNCTURE_ROOT/source-axisym/python:$PUNCTURE_ROOT/source-polar/examples:$PUNCTURE_ROOT/workflow-axisym/validation"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/horizon-v3-job-id.txt"
python3 "$PUNCTURE_ROOT/athenak-affine/tst/test_suite/z4c/check_hispid_controls.py" --executable "$PUNCTURE_ROOT/build-athenak-affine/src/athena" --manifest "$PUNCTURE_RUN/manifest.json" --cases gamma10 --affine-chart --boost-levels 4,8,12 --boost-flow-alpha .2 --flow-iterations 1000 --worker-timeout 300 --harmonic-storage factorized --output "$PUNCTURE_RUN/horizons-v3" > "$PUNCTURE_RUN/horizons-v3.log" 2>&1
