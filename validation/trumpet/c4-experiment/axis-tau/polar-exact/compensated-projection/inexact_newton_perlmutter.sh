#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/axis-tau-compensated-v3"
mkdir "$PUNCTURE_RUN"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
printf '%s  %s\n' 7398b4b254371c39409d511ebe0808fef5e8c45b9c1e6b1af2ccb52acecf31ea "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" | sha256sum -c -
export HISPID_TRACE_NEWTON=1
python3 "$PUNCTURE_ROOT/source-axis-tau-compensated/validation/run_trumpet_pilot.py" \
 --library "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" \
 --case moderate --n 240 --npolar 480 --nphi 16 --krylov-restart 64 \
 --memory-mib 65536 --linear-rtol .5 \
 --initial "$PUNCTURE_ROOT/axis-tau-compensated-v2/moderate240/diagnostic.checkpoint" \
 --output "$PUNCTURE_RUN/moderate240" > "$PUNCTURE_RUN/solve240.log" 2>&1
python3 - "$PUNCTURE_RUN/moderate240/result.json" <<'CHECK'
import json,sys
r=json.load(open(sys.argv[1]))
assert r['completed'] and r['diagnostics']['converged'],r['diagnostics']
assert r['minimum_psi']>0 and r['minimum_metric_eigenvalue']>0 and r['exterior_stencil_unmodified']
for region in ('near','bulk'):
 p=r['physical'][region]
 assert p['H_rms']<1e-6 and p['M_rms']<1e-6 and p['H_max']<1e-4 and p['M_max']<1e-4,(region,p)
CHECK
