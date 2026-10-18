#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-dual-endpoint-20261018
PUNCTURE_BUILD=$PUNCTURE_ROOT/build-by-capacity-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/by-matched-observers-v2
mkdir -p "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-$1.txt"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
# Reuse this allocation sequentially for one independent GPU correctness check.
# The native BY solve has not started; no concurrent srun step or second GPU.
if [[ "$1" == 224 && -f "$PUNCTURE_ROOT/focus-tau-build-v2/images.sha256" && ! -d "$PUNCTURE_ROOT/focus-tau-device-preflight-v1" ]]; then
 bash "$PUNCTURE_ROOT/focus_preflight_perlmutter.sh" || true
fi
python3 "$PUNCTURE_SOURCE/validation/compare_fixed_by.py" --library "$PUNCTURE_BUILD/libTwoPunctures.so" --hispid-result "$PUNCTURE_ROOT/qi-regularity-wide-v1/moderate$1/result.json" --observer "$PUNCTURE_ROOT/qi-regularity-wide-v1/moderate$1/physical.npz" --output "$PUNCTURE_RUN/by$1" > "$PUNCTURE_RUN/by$1.log" 2>&1
