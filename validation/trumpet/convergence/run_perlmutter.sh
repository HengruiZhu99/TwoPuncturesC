#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export PYTHONPATH="$PUNCTURE_ROOT/source-polar/python"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_ROOT/convergence/job-id.txt"
PUNCTURE_INITIAL="$PUNCTURE_ROOT/polar-refinement/moderate256x384/diagnostic.checkpoint"
for PUNCTURE_N in 192 224 256; do
  PUNCTURE_OUTPUT="$PUNCTURE_ROOT/convergence/moderate$PUNCTURE_N"
  python3 "$PUNCTURE_ROOT/workflow-polar/validation/run_trumpet_sequence_level.py" --library "$PUNCTURE_ROOT/build-polar-cuda/libHiSpID.so" --n "$PUNCTURE_N" --nphi 16 --krylov-restart 64 --memory-mib 65536 --initial "$PUNCTURE_INITIAL" --output "$PUNCTURE_OUTPUT" > "$PUNCTURE_ROOT/convergence/moderate$PUNCTURE_N.log" 2>&1
  python3 -c 'import json,sys;d=json.load(open(sys.argv[1]));print(sys.argv[1],d["physical"],flush=True);sys.exit(not(d["completed"] and d["diagnostics"]["converged"] and d["minimum_psi"]>0))' "$PUNCTURE_OUTPUT/result.json"
  PUNCTURE_INITIAL="$PUNCTURE_OUTPUT/diagnostic.checkpoint"
done
python3 "$PUNCTURE_ROOT/workflow-polar/validation/assess_trumpet_sequence.py" --runs "$PUNCTURE_ROOT/convergence/moderate192" "$PUNCTURE_ROOT/convergence/moderate224" "$PUNCTURE_ROOT/convergence/moderate256" --output "$PUNCTURE_ROOT/convergence/assessment.json"
