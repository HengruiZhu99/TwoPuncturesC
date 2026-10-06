#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/affine-horizon/binary256"
mkdir -p "$PUNCTURE_RUN"
export PYTHONPATH="$PUNCTURE_ROOT/source-axisym/python:$PUNCTURE_ROOT/workflow-axisym/validation:$PUNCTURE_ROOT/workflow-axisym/examples"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
python3 "$PUNCTURE_ROOT/workflow-axisym/validation/portable_sampler_migration.py" --checkpoint "$PUNCTURE_ROOT/axisymmetric/gamma10_256x320/diagnostic.checkpoint" --producer-library "$PUNCTURE_ROOT/build-axisym-cuda/libHiSpID.so" --consumer-library "$PUNCTURE_ROOT/build-focused-sampler/libHiSpID.so" --output "$PUNCTURE_RUN/migration.json" > "$PUNCTURE_RUN/migration.log" 2>&1
python3 "$PUNCTURE_ROOT/athenak-affine/tst/test_suite/z4c/check_hispid_binary.py" --executable "$PUNCTURE_ROOT/build-athenak-affine/src/athena" --checkpoint "$PUNCTURE_ROOT/axisymmetric/gamma10_256x320/diagnostic.checkpoint" --allow-diagnostic --migration-proof "$PUNCTURE_RUN/migration.json" --affine-chart --reuse-shapes --domain-half-width 16 --harmonic-storage factorized --geometry-threads 16 --strict-expansion --flow-alpha .2 --flow-iterations 1000 --levels 8,12,16 --timeout 1500 --output "$PUNCTURE_RUN/surfaces" > "$PUNCTURE_RUN/horizons.log" 2>&1
