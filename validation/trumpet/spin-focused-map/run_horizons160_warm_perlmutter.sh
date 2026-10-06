#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/spin-focused-map/horizon160-warm"
mkdir -p "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export PYTHONPATH="$PUNCTURE_ROOT/source-polar/python:$PUNCTURE_ROOT/workflow-polar/validation:$PUNCTURE_ROOT/workflow-polar/examples"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python3 "$PUNCTURE_ROOT/athenak-warm/tst/test_suite/z4c/check_hispid_binary.py" --executable "$PUNCTURE_ROOT/build-athenak-warm/src/athena" --checkpoint "$PUNCTURE_ROOT/spin-focused-map/spin99_160/diagnostic.checkpoint" --allow-diagnostic --migration-proof "$PUNCTURE_ROOT/spin-focused-map/horizon160/migration.json" --harmonic-storage factorized --geometry-threads 16 --strict-expansion --flow-alpha .2 --levels 12,16,20 --reuse-shapes --initial-shapes "$PUNCTURE_ROOT/spin-focused-map/horizon160-angular/surfaces/spectral_l12_n24/hispid.horizon_shape_0.txt" "$PUNCTURE_ROOT/spin-focused-map/horizon160-angular/surfaces/spectral_l12_n24/hispid.horizon_shape_1.txt" --timeout 1500 --output "$PUNCTURE_RUN/surfaces" > "$PUNCTURE_RUN/horizons.log" 2>&1
