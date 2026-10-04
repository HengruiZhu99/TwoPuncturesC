#!/usr/bin/env bash
set -euo pipefail
cd /pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/hs99uu-physical-d115c01
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=16 OMP_PROC_BIND=close OMP_PLACES=cores
PUNCTURE_PHYSICS_ROOT=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002
PUNCTURE_HS_OUTPUT="$PUNCTURE_PHYSICS_ROOT/hs99uu-physics-20261003"
scontrol show job "$SLURM_JOB_ID" > "$PUNCTURE_HS_OUTPUT/allocation.txt"
for PUNCTURE_HS_INDEX in 0 1 2; do
  PUNCTURE_HS_STATUS=0
  srun --mpi=none --resv-ports=0 -n1 -c32 --gpus=1 --cpu-bind=cores \
    python3 validation/run_extreme_kokkos.py \
    --library "$PUNCTURE_PHYSICS_ROOT/build-cuda/libHiSpID.so" \
    --performance-results "$PUNCTURE_PHYSICS_ROOT/human-override-performance-preserved-20261003.json" \
    --performance-artifact-root "$PUNCTURE_PHYSICS_ROOT/source" \
    --human-override validation/human_override_20261003.json \
    --plan validation/extreme_hs99uu_plan_20261003.json \
    --seed-controls "$PUNCTURE_PHYSICS_ROOT/extreme-human-override-20261003/seed-controls.json" \
    --source-floor-controls "$PUNCTURE_PHYSICS_ROOT/extreme-human-override-20261003/source-floor-spin99.json" \
    --case aligned_spin99_kokkos --grid-index "$PUNCTURE_HS_INDEX" \
    --threads 16 --geometry host --allow-diagnostic-investigation \
    --output-directory "$PUNCTURE_HS_OUTPUT/results" \
    > "$PUNCTURE_HS_OUTPUT/grid$PUNCTURE_HS_INDEX.log" 2>&1 || PUNCTURE_HS_STATUS=$?
  printf '%s\n' "$PUNCTURE_HS_STATUS" > "$PUNCTURE_HS_OUTPUT/grid$PUNCTURE_HS_INDEX.exit"
  python3 - "$PUNCTURE_HS_OUTPUT/results/results.json" "$PUNCTURE_HS_INDEX" <<'CHECK'
import json,sys
result=json.load(open(sys.argv[1]));index=int(sys.argv[2]);rows=result['aligned_spin99_kokkos_hs99uu']['records']
if len(rows)!=index+1:raise RuntimeError('new physical row did not complete; preserve logs and stop')
r=rows[-1]
print(json.dumps(dict(grid=r['resolution'],internally_converged=r['stopping_verified'],near_H=r['near']['H_rms'],near_M=r['near']['M_rms'],bulk_H=r['bulk']['H_rms'],bulk_M=r['bulk']['M_rms'],passed_strict=r['passed_strict'])),flush=True)
CHECK
 done
printf '%s\n' 'three_HS99UU_physical_rows_retained' > "$PUNCTURE_HS_OUTPUT/phase.txt"
