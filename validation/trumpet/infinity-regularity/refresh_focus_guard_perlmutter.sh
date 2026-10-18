#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-focus-tau-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/focus-tau-build-v2
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake --build "$PUNCTURE_ROOT/build-focus-tau-cuda" -j8 --target HiSpID test_hispid_focus test_hispid_kokkos > "$PUNCTURE_RUN/build-cuda.log" 2>&1
"$PUNCTURE_ROOT/build-focus-tau-cuda/test_hispid_focus" > "$PUNCTURE_RUN/focus-test.log"
cmake --build "$PUNCTURE_ROOT/build-focus-tau-sampler" -j8 --target HiSpID > "$PUNCTURE_RUN/build-sampler.log" 2>&1
python3 "$PUNCTURE_SOURCE/tests/test_parallel_observer.py" --library "$PUNCTURE_ROOT/build-focus-tau-sampler/libHiSpID.so" --threads 16 --output "$PUNCTURE_RUN/observer-parity.json" > "$PUNCTURE_RUN/observer-parity.log" 2>&1
cmake -S "$PUNCTURE_ROOT/athenak-infinity2-20261018" -B "$PUNCTURE_ROOT/build-athenak-focus-tau" > "$PUNCTURE_RUN/configure-athena.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-athenak-focus-tau" -j16 --target athena > "$PUNCTURE_RUN/build-athena.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-focus-tau-cuda/"lib*.so "$PUNCTURE_ROOT/build-focus-tau-sampler/"lib*.so "$PUNCTURE_ROOT/build-athenak-focus-tau/src/athena" > "$PUNCTURE_RUN/images.sha256"
