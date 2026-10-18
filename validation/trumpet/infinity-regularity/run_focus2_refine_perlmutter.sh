#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-focus2-20261018
PUNCTURE_BUILD=$PUNCTURE_ROOT/build-focus2-cuda
PUNCTURE_RUN=$PUNCTURE_ROOT/focus2-refine-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores HISPID_TRACE_NEWTON=1
sha256sum -c "$PUNCTURE_ROOT/focus2-build-v3/images.sha256" > "$PUNCTURE_RUN/preflight-images.log"
cp "$PUNCTURE_ROOT/focus2-v2/actual-kernels.log" "$PUNCTURE_RUN/actual-kernels.log"
sha256sum "$PUNCTURE_BUILD/"lib*.so > "$PUNCTURE_RUN/images.sha256"
PUNCTURE_INITIAL=()
for PUNCTURE_N in 192 224 256; do
 if [[ "$PUNCTURE_N" == 192 ]]; then
 PUNCTURE_INITIAL=(--initial "$PUNCTURE_ROOT/focus-tau-v1/moderate192/diagnostic.checkpoint" --initial-source-library "$PUNCTURE_ROOT/build-focus-dd-cuda/libHiSpID.so")
 fi
 python3 "$PUNCTURE_SOURCE/validation/run_trumpet_pilot.py" \
  --library "$PUNCTURE_BUILD/libHiSpID.so" --case qi-moderate-wide \
  --n "$PUNCTURE_N" --npolar "$((2*PUNCTURE_N))" --nphi 16 \
  --krylov-restart 64 --max-krylov 512 --linear-rtol .1 \
  --observer-threads 16 --memory-mib 65536 --skip-observer-on-failure "${PUNCTURE_INITIAL[@]}" \
  --output "$PUNCTURE_RUN/moderate$PUNCTURE_N" > "$PUNCTURE_RUN/solve$PUNCTURE_N.log" 2>&1
 python3 - "$PUNCTURE_RUN/moderate$PUNCTURE_N/result.json" <<'CHECK'
import json,sys
r=json.load(open(sys.argv[1]));assert r['completed'] and r['diagnostics']['converged'],r['diagnostics']
assert r['minimum_psi']>0 and r['minimum_metric_eigenvalue']>0 and r['exterior_stencil_unmodified']
for region in ('near','bulk'):
 q=r['physical'][region];assert q['H_rms']<1e-6 and q['M_rms']<1e-6 and q['H_max']<1e-4 and q['M_max']<1e-4,(region,q)
if r['config']['n'][0]==192:
 assert r['physical']['bulk']['M_rms']<1.0538120185710488e-6/4 and r['physical']['bulk']['H_rms']<1.9591214219981626e-9,'no clear refinement improvement'
CHECK
 PUNCTURE_INITIAL=(--initial "$PUNCTURE_RUN/moderate$PUNCTURE_N/diagnostic.checkpoint")
done
PUNCTURE_SEQUENCE_STATUS=0
python3 "$PUNCTURE_SOURCE/validation/assess_trumpet_sequence.py" \
 --runs "$PUNCTURE_RUN/moderate192" "$PUNCTURE_RUN/moderate224" "$PUNCTURE_RUN/moderate256" \
 --output "$PUNCTURE_RUN/assessment.json" > "$PUNCTURE_RUN/assessment.log" 2>&1 || PUNCTURE_SEQUENCE_STATUS=$?
exit "$PUNCTURE_SEQUENCE_STATUS"
