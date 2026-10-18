#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_CASE=${1:-qi-moderate}
case "$PUNCTURE_CASE" in
 qi-moderate) PUNCTURE_SOURCE="$PUNCTURE_ROOT/workflow-regularity-control-20261018"; PUNCTURE_RUN="$PUNCTURE_ROOT/qi-regularity-control-v1";;
 qi-moderate-wide) PUNCTURE_SOURCE="$PUNCTURE_ROOT/workflow-regularity-wide-20261018"; PUNCTURE_RUN="$PUNCTURE_ROOT/qi-regularity-wide-v1";;
 *) exit 2;;
esac
PUNCTURE_LIBRARY="$PUNCTURE_ROOT/build-arnoldi-20261018/libHiSpID.so"
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
printf '%s  %s\n' c4c45eaee847a6ab17ad794c0d81321ea27d2a2cc313f84c8f0d622e7d8b4741 "$PUNCTURE_LIBRARY" | sha256sum -c -
printf '%s  %s\n' e4ecc26cbffb4ac2e18574ac5c6e7e5b6086f49fee8071e69d3eac90ed03db8e "$PUNCTURE_ROOT/build-arnoldi-20261018/libTwoPunctures.so" | sha256sum -c -
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores HISPID_TRACE_NEWTON=1
PUNCTURE_INITIAL=()
for PUNCTURE_N in 160 192 224; do
 python3 "$PUNCTURE_SOURCE/validation/run_trumpet_pilot.py" \
  --library "$PUNCTURE_LIBRARY" --case "$PUNCTURE_CASE" \
  --n "$PUNCTURE_N" --npolar "$((2*PUNCTURE_N))" --nphi 16 \
  --krylov-restart 64 --max-krylov 512 --linear-rtol .1 \
  --memory-mib 65536 --skip-observer-on-failure "${PUNCTURE_INITIAL[@]}" \
  --output "$PUNCTURE_RUN/moderate$PUNCTURE_N" > "$PUNCTURE_RUN/solve$PUNCTURE_N.log" 2>&1
 python3 - "$PUNCTURE_RUN/moderate$PUNCTURE_N/result.json" <<'CHECK'
import json,sys
r=json.load(open(sys.argv[1]))
assert r['completed'] and r['diagnostics']['converged'],r['diagnostics']
assert r['minimum_psi']>0 and r['minimum_metric_eigenvalue']>0 and r['exterior_stencil_unmodified']
for region in ('near','bulk'):
 q=r['physical'][region]
 assert q['H_rms']<1e-6 and q['M_rms']<1e-6 and q['H_max']<1e-4 and q['M_max']<1e-4,(region,q)
CHECK
 PUNCTURE_INITIAL=(--initial "$PUNCTURE_RUN/moderate$PUNCTURE_N/diagnostic.checkpoint")
done
python3 "$PUNCTURE_SOURCE/validation/assess_trumpet_sequence.py" \
 --runs "$PUNCTURE_RUN/moderate160" "$PUNCTURE_RUN/moderate192" "$PUNCTURE_RUN/moderate224" \
 --output "$PUNCTURE_RUN/assessment.json" > "$PUNCTURE_RUN/assessment.log" 2>&1
