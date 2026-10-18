#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN=$PUNCTURE_ROOT/focus-tau-device-preflight-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
PUNCTURE_BUILD=$PUNCTURE_ROOT/build-focus-tau-cuda
sha256sum "$PUNCTURE_BUILD/test_hispid_kokkos" "$PUNCTURE_BUILD/"lib*.so > "$PUNCTURE_RUN/images.sha256"
"$PUNCTURE_BUILD/test_hispid_kokkos" > "$PUNCTURE_RUN/actual-kernels.log" 2>&1
touch "$PUNCTURE_RUN/passed"
