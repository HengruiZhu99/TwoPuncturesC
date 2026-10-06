#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/axis-tau-polar-sequence-v1"
mkdir "$PUNCTURE_RUN"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores HISPID_TRACE_NEWTON=1
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
printf '%s  %s\n' 9489bca44c0bbd6ea20ae4996b4e2cdf225e84e3436ec80b253c4ffaf2b0c3e6 "$PUNCTURE_ROOT/build-axis-tau-polar/libHiSpID.so" | sha256sum -c -
PUNCTURE_INITIAL="$PUNCTURE_ROOT/axis-tau-polar-full-v1/moderate192/diagnostic.checkpoint"
for PUNCTURE_N in 224 240 256; do
 PUNCTURE_NPOLAR=$((2*PUNCTURE_N))
 python3 "$PUNCTURE_ROOT/source-axis-tau-polar/validation/run_trumpet_pilot.py" \
  --library "$PUNCTURE_ROOT/build-axis-tau-polar/libHiSpID.so" \
  --case moderate --n "$PUNCTURE_N" --npolar "$PUNCTURE_NPOLAR" --nphi 16 \
  --krylov-restart 24 --memory-mib 65536 --linear-rtol .001 \
  --initial "$PUNCTURE_INITIAL" --output "$PUNCTURE_RUN/moderate$PUNCTURE_N" > "$PUNCTURE_RUN/solve$PUNCTURE_N.log" 2>&1
 python3 - "$PUNCTURE_RUN/moderate$PUNCTURE_N/result.json" <<'CHECK'
import json,sys
x=json.load(open(sys.argv[1]))
assert x['completed'] and x['diagnostics']['converged']
assert x['minimum_psi']>0 and x['minimum_metric_eigenvalue']>0 and x['exterior_stencil_unmodified']
for region in ('near','bulk'):
 p=x['physical'][region]
 assert p['H_rms']<1e-6 and p['M_rms']<1e-6 and p['H_max']<1e-4 and p['M_max']<1e-4,(region,p)
print(sys.argv[1],x['physical'],flush=True)
CHECK
 PUNCTURE_INITIAL="$PUNCTURE_RUN/moderate$PUNCTURE_N/diagnostic.checkpoint"
done
python3 "$PUNCTURE_ROOT/source-axis-tau-polar/validation/assess_trumpet_sequence.py" \
 --runs "$PUNCTURE_RUN/moderate224" "$PUNCTURE_RUN/moderate240" "$PUNCTURE_RUN/moderate256" \
 --output "$PUNCTURE_RUN/assessment.json" > "$PUNCTURE_RUN/assessment.log" 2>&1
