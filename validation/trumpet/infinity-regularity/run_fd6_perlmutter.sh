#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-higher-observer-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/focus-fd6-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python3 "$PUNCTURE_SOURCE/tests/test_physical_fd6.py" --points "$PUNCTURE_ROOT/focus-tau-v1/moderate160/physical.npz" --output "$PUNCTURE_RUN/control.json" > "$PUNCTURE_RUN/control.log" 2>&1
for PUNCTURE_N in 160 192 224; do
 for PUNCTURE_WAIT in $(seq 1 60); do
  if python3 - "$PUNCTURE_ROOT/focus-tau-v1/moderate$PUNCTURE_N/result.json" <<'READY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);raise SystemExit(0 if p.exists() and json.loads(p.read_text()).get('completed') else 1)
READY
  then break; fi
  sleep 10
 done
 python3 "$PUNCTURE_SOURCE/validation/observe_saved_fd6.py" --system hispid --run "$PUNCTURE_ROOT/focus-tau-v1/moderate$PUNCTURE_N" --library "$PUNCTURE_ROOT/build-focus-dd-cuda/libHiSpID.so" --output "$PUNCTURE_RUN/hispid$PUNCTURE_N" > "$PUNCTURE_RUN/hispid$PUNCTURE_N.log" 2>&1
 python3 "$PUNCTURE_SOURCE/validation/observe_saved_fd6.py" --system by --run "$PUNCTURE_ROOT/by-matched-observers-v2/by$PUNCTURE_N" --library "$PUNCTURE_ROOT/build-by-capacity-20261018/libTwoPunctures.so" --output "$PUNCTURE_RUN/by$PUNCTURE_N" > "$PUNCTURE_RUN/by$PUNCTURE_N.log" 2>&1
done
python3 "$PUNCTURE_SOURCE/validation/assess_saved_fd6.py" --runs "$PUNCTURE_RUN" --source "$PUNCTURE_ROOT/focus-tau-v1" --output "$PUNCTURE_RUN/assessment.json" > "$PUNCTURE_RUN/assessment.log" 2>&1
