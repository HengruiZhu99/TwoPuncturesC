#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE="$PUNCTURE_ROOT/source-modal-relaxation-20261018"
PUNCTURE_BUILD="$PUNCTURE_ROOT/build-modal-relaxation-20261018"
PUNCTURE_RUN="$PUNCTURE_ROOT/modal-relaxation-sequence-v1"
PUNCTURE_MIDDLE="$PUNCTURE_ROOT/modal-relaxation-pilot-v1/moderate240"
python3 - "$PUNCTURE_MIDDLE/result.json" <<'CHECK'
import json,sys
r=json.load(open(sys.argv[1]));assert r['completed'] and r['diagnostics']['converged']
assert r['exterior_stencil_unmodified'] and r['minimum_psi']>0 and r['minimum_metric_eigenvalue']>0
for region in ('near','bulk'):
 q=r['physical'][region];assert q['H_rms']<1e-6 and q['M_rms']<1e-6 and q['H_max']<1e-4 and q['M_max']<1e-4
CHECK
printf '%s  %s\n' 7abe943bb14d35b49dc79599431abf4cd7666eb072ed317387009e49c3fb9c03 "$PUNCTURE_BUILD/libHiSpID.so" | sha256sum -c -
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores HISPID_TRACE_NEWTON=1
for PUNCTURE_N in 224 256; do
 if [ "$PUNCTURE_N" = 224 ]; then
  PUNCTURE_INITIAL="$PUNCTURE_ROOT/axis-tau-polar-sequence-v2/moderate224/diagnostic.checkpoint"
  PUNCTURE_EXTRA=(--initial-source-library "$PUNCTURE_ROOT/build-axis-tau-polar/libHiSpID.so")
 else
  PUNCTURE_INITIAL="$PUNCTURE_MIDDLE/diagnostic.checkpoint"
  PUNCTURE_EXTRA=()
 fi
 python3 "$PUNCTURE_SOURCE/validation/run_trumpet_pilot.py" \
  --library "$PUNCTURE_BUILD/libHiSpID.so" --case moderate --n "$PUNCTURE_N" \
  --npolar "$((2*PUNCTURE_N))" --nphi 16 --krylov-restart 64 --max-krylov 512 \
  --memory-mib 65536 --linear-rtol .1 --initial "$PUNCTURE_INITIAL" "${PUNCTURE_EXTRA[@]}" \
  --output "$PUNCTURE_RUN/moderate$PUNCTURE_N" > "$PUNCTURE_RUN/solve$PUNCTURE_N.log" 2>&1
 python3 - "$PUNCTURE_RUN/moderate$PUNCTURE_N/result.json" <<'CHECK'
import json,sys
r=json.load(open(sys.argv[1]));assert r['completed'] and r['diagnostics']['converged'],r['diagnostics']
assert r['minimum_psi']>0 and r['minimum_metric_eigenvalue']>0 and r['exterior_stencil_unmodified']
for region in ('near','bulk'):
 q=r['physical'][region];assert q['H_rms']<1e-6 and q['M_rms']<1e-6 and q['H_max']<1e-4 and q['M_max']<1e-4,(region,q)
CHECK
done
python3 "$PUNCTURE_SOURCE/validation/assess_trumpet_sequence.py" \
 --runs "$PUNCTURE_RUN/moderate224" "$PUNCTURE_MIDDLE" "$PUNCTURE_RUN/moderate256" \
 --output "$PUNCTURE_RUN/assessment.json" > "$PUNCTURE_RUN/assessment.log" 2>&1
