#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN=$PUNCTURE_ROOT/focus-dd-build-v3
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake --build "$PUNCTURE_ROOT/build-focus-dd-cuda" -j8 --target HiSpID test_hispid_focus test_hispid_kokkos > "$PUNCTURE_RUN/build.log" 2>&1
"$PUNCTURE_ROOT/build-focus-dd-cuda/test_hispid_focus" > "$PUNCTURE_RUN/focus-test.log"
sha256sum "$PUNCTURE_ROOT/build-focus-dd-cuda/"lib*.so "$PUNCTURE_ROOT/build-focus-dd-cuda/test_hispid_kokkos" > "$PUNCTURE_RUN/images.sha256"
