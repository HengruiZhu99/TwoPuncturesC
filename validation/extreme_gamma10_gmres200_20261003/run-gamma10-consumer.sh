#!/usr/bin/env bash
set -euo pipefail
PUNCTURE_PHYSICS_ROOT=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002
PUNCTURE_GAMMA_OUTPUT="$PUNCTURE_PHYSICS_ROOT/gamma10-gmres200-physics-20261003"
PUNCTURE_GAMMA_CHECKPOINT="$PUNCTURE_GAMMA_OUTPUT/results/headon_gamma10_kokkos_gmres200_128_256_8.checkpoint"
cd "$PUNCTURE_PHYSICS_ROOT/gamma10-physical-cacdace"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores
srun --mpi=none --resv-ports=0 -n1 -c1 --gpus=1 --cpu-bind=cores \
  python3 validation/portable_sampler_migration.py \
  --checkpoint "$PUNCTURE_GAMMA_CHECKPOINT" \
  --producer-library "$PUNCTURE_PHYSICS_ROOT/build-cuda/libHiSpID.so" \
  --consumer-library "$PUNCTURE_PHYSICS_ROOT/build-reference/libHiSpID.so" \
  --timeout 600 --output "$PUNCTURE_GAMMA_OUTPUT/consumer-proof/gamma10-gmres200.json" \
  > "$PUNCTURE_GAMMA_OUTPUT/consumer-proof.log" 2>&1
export PYTHONPATH="$PUNCTURE_PHYSICS_ROOT/athenak-030cbcc/tst/test_suite/z4c:$PUNCTURE_PHYSICS_ROOT/athenak-import-f231d722"
srun --mpi=none --resv-ports=0 -n1 -c1 --gpus=1 --cpu-bind=cores \
  python3 "$PUNCTURE_PHYSICS_ROOT/athenak-import-f231d722/check_hispid_import.py" \
  --executable "$PUNCTURE_PHYSICS_ROOT/athenak-build-extreme-serial-20261003/src/athena" \
  --template "$PUNCTURE_PHYSICS_ROOT/athenak-030cbcc/tst/inputs/hispid.athinput" \
  --checkpoint "$PUNCTURE_GAMMA_CHECKPOINT" --allow-diagnostic --domain-half-width 40 --timeout 600 \
  --migration-proof "$PUNCTURE_GAMMA_OUTPUT/consumer-proof/gamma10-gmres200.json" \
  --output "$PUNCTURE_GAMMA_OUTPUT/import-only" \
  > "$PUNCTURE_GAMMA_OUTPUT/import-only.log" 2>&1
printf '%s\n' 'Gamma10_consumer_checks_terminal' > "$PUNCTURE_GAMMA_OUTPUT/consumer-phase.txt"
