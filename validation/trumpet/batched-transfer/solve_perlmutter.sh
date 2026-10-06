#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
# Contingency only: inspect the terminal result of job59413586 before submitting.
# If it produced usable radial512 data, use that data instead of rerunning it.
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/batched-transfer-gamma512-v1"
PUNCTURE_LIBRARY="$PUNCTURE_ROOT/build-batched-transfer-cuda/libHiSpID.so"
printf '%s  %s\n' fdaf287805d5f66dc66a72cd9b445188bfa1a1834d818c4eacfc270fedbd1b61 "$PUNCTURE_LIBRARY" | sha256sum --check
# A fresh directory preserves every earlier run, including any timeout evidence.
mkdir "$PUNCTURE_RUN"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
export PYTHONPATH="$PUNCTURE_ROOT/source-batched-transfer/python:$PUNCTURE_ROOT/source-batched-transfer/validation:$PUNCTURE_ROOT/source-batched-transfer/examples"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
python3 "$PUNCTURE_ROOT/source-batched-transfer/validation/run_trumpet_pilot.py" \
  --library "$PUNCTURE_LIBRARY" --case gamma10 --axisymmetric \
  --n 512 --npolar 320 --nphi 8 --krylov-restart 64 --memory-mib 65536 \
  --initial "$PUNCTURE_ROOT/axisymmetric/gamma10_256x320/diagnostic.checkpoint" \
  --initial-source-library "$PUNCTURE_ROOT/build-axisym-cuda/libHiSpID.so" \
  --output "$PUNCTURE_RUN/solve" > "$PUNCTURE_RUN/solve.log" 2>&1
