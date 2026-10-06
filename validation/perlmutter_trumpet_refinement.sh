#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?allocated compute node required}"
PUNCTURE_ROOT=${1:?run root}
export PYTHONPATH="$PUNCTURE_ROOT/source/python:$PUNCTURE_ROOT/workflow-refinement/examples:$PUNCTURE_ROOT/source/validation"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=16 OMP_PROC_BIND=spread OMP_PLACES=cores
PUNCTURE_OUT="$PUNCTURE_ROOT/refinement-physics"
mkdir "$PUNCTURE_OUT"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_OUT/job-id.txt"
python3 "$PUNCTURE_ROOT/source/validation/check_trumpet_near_seed_refined.py" --library "$PUNCTURE_ROOT/build-sampler/libHiSpID.so" --output "$PUNCTURE_OUT/near-spin99-refined.json" --case spin99 --factors .005,.0025 > "$PUNCTURE_OUT/near-spin99.log" 2>&1
python3 "$PUNCTURE_ROOT/workflow-refinement/validation/inspect_trumpet_checkpoint.py" --library "$PUNCTURE_ROOT/build-cuda/libHiSpID.so" --pilot "$PUNCTURE_ROOT/first-physics/moderate-pilot" --output "$PUNCTURE_OUT/coarse-observer.json" > "$PUNCTURE_OUT/coarse-observer.log" 2>&1
python3 "$PUNCTURE_ROOT/workflow-refinement/validation/run_trumpet_pilot.py" --library "$PUNCTURE_ROOT/build-cuda/libHiSpID.so" --output "$PUNCTURE_OUT/moderate40" --n 40 --nphi 12 --initial "$PUNCTURE_ROOT/first-physics/moderate-pilot/diagnostic.checkpoint" > "$PUNCTURE_OUT/moderate40.log" 2>&1
printf 'completed\n' > "$PUNCTURE_OUT/completed.txt"
