#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/axial-factors/consumer-v1"
# Run only after reviewing the completed radial-512 physical diagnostics.
# The driver verifies the checkpoint-bound sampler proof and stops on a failed row.
# No failed radial-256 surface is reused as a qualified measurement.
mkdir "$PUNCTURE_RUN"
export PYTHONPATH="$PUNCTURE_ROOT/source-axial-factors/python:$PUNCTURE_ROOT/workflow-radial512/validation:$PUNCTURE_ROOT/workflow-radial512/examples"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
python3 "$PUNCTURE_ROOT/workflow-radial512/validation/portable_sampler_migration.py" --checkpoint "$PUNCTURE_ROOT/axial-factors/gamma10_512x320-v1/diagnostic.checkpoint" --producer-library "$PUNCTURE_ROOT/build-axial-factors-cuda/libHiSpID.so" --consumer-library "$PUNCTURE_ROOT/build-radial512-sampler/libHiSpID.so" --output "$PUNCTURE_RUN/migration.json" > "$PUNCTURE_RUN/migration.log" 2>&1
python3 "$PUNCTURE_ROOT/athenak-affine/tst/test_suite/z4c/check_hispid_binary.py" --executable "$PUNCTURE_ROOT/build-athenak-radial512/src/athena" --checkpoint "$PUNCTURE_ROOT/axial-factors/gamma10_512x320-v1/diagnostic.checkpoint" --allow-diagnostic --migration-proof "$PUNCTURE_RUN/migration.json" --affine-chart --reuse-shapes --domain-half-width 16 --harmonic-storage factorized --geometry-threads 16 --strict-expansion --flow-alpha .2 --flow-iterations 1000 --levels 8,12,16 --timeout 1500 --output "$PUNCTURE_RUN/surfaces" > "$PUNCTURE_RUN/horizons.log" 2>&1
