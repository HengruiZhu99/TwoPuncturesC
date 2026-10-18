#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE="$PUNCTURE_ROOT/workflow-moderate-wide-20261018"
PUNCTURE_RUN="$PUNCTURE_ROOT/moderate-wide-fine-v1"
PUNCTURE_COARSE="$PUNCTURE_ROOT/moderate-wide-sequence-v1/moderate192"
mkdir "$PUNCTURE_RUN"
printf '%s  %s\n' 40d6a2165f07298427aa8acbca4d808748888eb8a579febe15815d02622f574c "$PUNCTURE_ROOT/build-axis-tau-compensated/libTwoPunctures.so" | sha256sum -c -
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores HISPID_TRACE_NEWTON=1
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
printf '%s  %s\n' 7398b4b254371c39409d511ebe0808fef5e8c45b9c1e6b1af2ccb52acecf31ea "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" | sha256sum -c -
python3 - "$PUNCTURE_COARSE/result.json" <<'CHECK'
import json,sys
r=json.load(open(sys.argv[1]));assert r['completed'] and r['diagnostics']['converged'] and r['exterior_stencil_unmodified']
for region in ('near','bulk'):
 q=r['physical'][region]
 assert q['H_rms']<1e-6 and q['M_rms']<1e-6 and q['H_max']<1e-4 and q['M_max']<1e-4,(region,q)
CHECK
PUNCTURE_INITIAL="$PUNCTURE_ROOT/moderate-wide-sequence-v1/moderate224/diagnostic.checkpoint"
for PUNCTURE_N in 240; do
 python3 "$PUNCTURE_SOURCE/validation/run_trumpet_pilot.py" \
  --library "$PUNCTURE_ROOT/build-axis-tau-compensated/libHiSpID.so" \
  --case moderate-wide --n "$PUNCTURE_N" --npolar "$((2*PUNCTURE_N))" --nphi 16 \
  --krylov-restart 64 --memory-mib 65536 --linear-rtol .1 \
  --initial "$PUNCTURE_INITIAL" --output "$PUNCTURE_RUN/moderate$PUNCTURE_N" > "$PUNCTURE_RUN/solve$PUNCTURE_N.log" 2>&1
 python3 - "$PUNCTURE_RUN/moderate$PUNCTURE_N/result.json" <<'CHECK'
import json,sys
x=json.load(open(sys.argv[1]));assert x['completed'] and x['diagnostics']['converged'],x['diagnostics']
CHECK
 PUNCTURE_INITIAL="$PUNCTURE_RUN/moderate$PUNCTURE_N/diagnostic.checkpoint"
done
python3 "$PUNCTURE_SOURCE/validation/assess_trumpet_sequence.py" \
 --runs "$PUNCTURE_COARSE" "$PUNCTURE_ROOT/moderate-wide-sequence-v1/moderate224" "$PUNCTURE_RUN/moderate240" \
 --output "$PUNCTURE_RUN/assessment.json" > "$PUNCTURE_RUN/assessment.log" 2>&1
