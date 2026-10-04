#!/usr/bin/env bash
set -euo pipefail
cd /pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/gamma10-physical-cacdace
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=16 OMP_PROC_BIND=close OMP_PLACES=cores
PUNCTURE_PHYSICS_ROOT=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002
PUNCTURE_GAMMA_OUTPUT="$PUNCTURE_PHYSICS_ROOT/gamma10-gmres200-physics-20261003"
scontrol show job "$SLURM_JOB_ID" > "$PUNCTURE_GAMMA_OUTPUT/allocation.txt"
PUNCTURE_GAMMA_STATUS=0
srun --mpi=none --resv-ports=0 -n1 -c32 --gpus=1 --cpu-bind=cores \
  python3 validation/run_extreme_kokkos.py \
  --library "$PUNCTURE_PHYSICS_ROOT/build-cuda/libHiSpID.so" \
  --performance-results "$PUNCTURE_PHYSICS_ROOT/human-override-performance-preserved-20261003.json" \
  --performance-artifact-root "$PUNCTURE_PHYSICS_ROOT/source" \
  --human-override validation/human_override_20261003.json \
  --plan validation/extreme_gamma10_gmres200_plan_20261003.json \
  --seed-controls "$PUNCTURE_PHYSICS_ROOT/extreme-human-override-20261003/seed-controls.json" \
  --source-floor-controls "$PUNCTURE_PHYSICS_ROOT/extreme-human-override-20261003/source-floor-gamma10.json" \
  --initial-checkpoint "$PUNCTURE_PHYSICS_ROOT/extreme-human-override-20261003/results/headon_gamma10_kokkos_128_256_8.checkpoint" \
  --case headon_gamma10_kokkos --grid-index 0 --threads 16 --geometry host \
  --allow-diagnostic-investigation --output-directory "$PUNCTURE_GAMMA_OUTPUT/results" \
  > "$PUNCTURE_GAMMA_OUTPUT/trial.log" 2>&1 || PUNCTURE_GAMMA_STATUS=$?
printf '%s\n' "$PUNCTURE_GAMMA_STATUS" > "$PUNCTURE_GAMMA_OUTPUT/trial.exit"
printf '%s\n' 'Gamma10_trial_terminal' > "$PUNCTURE_GAMMA_OUTPUT/phase.txt"
