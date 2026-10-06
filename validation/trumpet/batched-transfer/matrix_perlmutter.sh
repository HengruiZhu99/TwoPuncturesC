#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/opt/cray/pe/gcc-native/14/bin/g++ -O2 -std=c++17 "$PUNCTURE_ROOT/batched-transfer-source/tests/test_hispid_modal_transfer.cpp" -lgsl -lgslcblas -o "$PUNCTURE_ROOT/batched-transfer/matrix-control" > "$PUNCTURE_ROOT/batched-transfer/build.log" 2>&1
"$PUNCTURE_ROOT/batched-transfer/matrix-control" > "$PUNCTURE_ROOT/batched-transfer/matrix.txt"
sha256sum "$PUNCTURE_ROOT/batched-transfer-source/tests/test_hispid_modal_transfer.cpp" "$PUNCTURE_ROOT/batched-transfer-source/src/HiSpID_modal_transfer.hpp" > "$PUNCTURE_ROOT/batched-transfer/source-manifest.sha256"
