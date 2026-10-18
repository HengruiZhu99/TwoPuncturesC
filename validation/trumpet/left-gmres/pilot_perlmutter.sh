#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE="$PUNCTURE_ROOT/source-left-gmres-20261018"
PUNCTURE_BUILD="$PUNCTURE_ROOT/build-left-gmres-20261018"
PUNCTURE_RUN="$PUNCTURE_ROOT/left-gmres-pilot-v1"
bash "$PUNCTURE_ROOT/left-gmres-build.sh"
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python3 "$PUNCTURE_SOURCE/validation/trumpet/c4-experiment/scaled-inverse/control.py" \
 --library "$PUNCTURE_BUILD/libHiSpID.so" --execution kokkos --solve \
 --output "$PUNCTURE_RUN/control.json" > "$PUNCTURE_RUN/control.log" 2>&1
python3 - "$PUNCTURE_RUN/control.json" "$PUNCTURE_SOURCE/validation/trumpet/left-gmres/control.json" <<'CHECK'
import json,sys,numpy as np
from pathlib import Path
r=json.load(open(sys.argv[1]));c=json.load(open(sys.argv[2]));assert r['diagnostics']['converged'] and c['diagnostics']['converged']
e=max(float(np.max(abs(np.array(r['field_witness'][k])-np.array(c['field_witness'][k]))/(1+abs(np.array(c['field_witness'][k]))))) for k in c['field_witness'])
Path(sys.argv[1]).with_name('backend-comparison.json').write_text(json.dumps(dict(field_scaled_difference=e,limit=1e-10))+'\n');assert e<1e-10
assert r['resolved_options']['preconditioner']=='modal_left_original_norm_confirmed'
CHECK
export HISPID_TRACE_NEWTON=1 PUNCTURE_KRYLOV_ARNOLDI="$PUNCTURE_RUN/arnoldi.jsonl"
python3 "$PUNCTURE_SOURCE/validation/run_trumpet_pilot.py" \
 --library "$PUNCTURE_BUILD/libHiSpID.so" --case moderate \
 --n 240 --npolar 480 --nphi 16 --krylov-restart 64 --max-krylov 512 \
 --memory-mib 65536 --linear-rtol .1 --skip-observer-on-failure \
 --initial "$PUNCTURE_ROOT/axis-tau-compensated-v3/moderate240/diagnostic.checkpoint" \
 --initial-source-library "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" \
 --output "$PUNCTURE_RUN/moderate240" > "$PUNCTURE_RUN/solve240.log" 2>&1
