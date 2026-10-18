#!/usr/bin/env bash
set -uo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-focus2-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/focus2-fd6-v1
mkdir "$PUNCTURE_RUN" || exit 1
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
cp "$PUNCTURE_ROOT/focus-fd6-v1/control.json" "$PUNCTURE_RUN/control.json" || exit 1
for PUNCTURE_N in 160 192 224; do
 python3 - "$PUNCTURE_ROOT/focus-tau-v1/moderate$PUNCTURE_N/physical.npz" "$PUNCTURE_ROOT/focus2-v2/moderate$PUNCTURE_N/physical.npz" <<'CHECK'
import sys,numpy as np
with np.load(sys.argv[1]) as old,np.load(sys.argv[2]) as new:
 assert np.array_equal(old['xyz'],new['xyz']) and np.array_equal(old['step'],new['step'])
CHECK
 if [[ "$?" != 0 ]]; then exit 1;fi
 python3 "$PUNCTURE_SOURCE/validation/observe_saved_fd6.py" --system hispid --run "$PUNCTURE_ROOT/focus2-v2/moderate$PUNCTURE_N" --library "$PUNCTURE_ROOT/build-focus2-cuda/libHiSpID.so" --output "$PUNCTURE_RUN/hispid$PUNCTURE_N" > "$PUNCTURE_RUN/hispid$PUNCTURE_N.log" 2>&1
 if [[ "$?" != 0 ]]; then exit 1;fi
done
python3 "$PUNCTURE_SOURCE/validation/assess_saved_fd6.py" --runs "$PUNCTURE_RUN" --source "$PUNCTURE_ROOT/focus2-v2" --output "$PUNCTURE_RUN/assessment.json" > "$PUNCTURE_RUN/assessment.log" 2>&1
PUNCTURE_FD_STATUS=$?
bash "$PUNCTURE_ROOT/horizons_focus2_perlmutter.sh"
PUNCTURE_HORIZON_STATUS=$?
printf '{"fd6_exit_code":%d,"horizon_exit_code":%d}\n' "$PUNCTURE_FD_STATUS" "$PUNCTURE_HORIZON_STATUS" > "$PUNCTURE_ROOT/focus2-completion-status.json"
if [[ "$PUNCTURE_FD_STATUS" != 0 || "$PUNCTURE_HORIZON_STATUS" != 0 ]]; then exit 1;fi
