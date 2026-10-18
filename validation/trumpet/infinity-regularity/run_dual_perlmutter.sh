#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-dual-endpoint-20261018
PUNCTURE_BUILD=$PUNCTURE_ROOT/build-dual-endpoint-cuda
PUNCTURE_RUN=$PUNCTURE_ROOT/dual-endpoint-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores HISPID_TRACE_NEWTON=1
python3 "$PUNCTURE_SOURCE/tests/check_geometry_backend.py" --library "$PUNCTURE_BUILD/libHiSpID.so" --output "$PUNCTURE_RUN/backend-test.json" > "$PUNCTURE_RUN/backend-test.log" 2>&1
sha256sum "$PUNCTURE_BUILD/"lib*.so > "$PUNCTURE_RUN/images.sha256"
PUNCTURE_INITIAL=()
for PUNCTURE_N in 160 192 224; do
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
CHECK
 if [ "$PUNCTURE_N" -gt 160 ]; then
  python3 - "$PUNCTURE_RUN" "$PUNCTURE_N" <<'PREFIX'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);n=int(sys.argv[2]);old=json.loads((p/f'moderate{n-32}/result.json').read_text());new=json.loads((p/f'moderate{n}/result.json').read_text())
checks={f'{r}_{k}':new['physical'][r][k]<old['physical'][r][k] for r in ('near','bulk') for k in ('H_rms','M_rms')}
(p/f'prefix{n}.json').write_text(json.dumps(dict(passed=all(checks.values()),checks=checks),indent=2)+'\n')
assert all(checks.values()),checks
PREFIX
 fi
 PUNCTURE_INITIAL=(--initial "$PUNCTURE_RUN/moderate$PUNCTURE_N/diagnostic.checkpoint")
done
python3 "$PUNCTURE_SOURCE/validation/assess_trumpet_sequence.py" \
 --runs "$PUNCTURE_RUN/moderate160" "$PUNCTURE_RUN/moderate192" "$PUNCTURE_RUN/moderate224" \
 --output "$PUNCTURE_RUN/assessment.json" > "$PUNCTURE_RUN/assessment.log" 2>&1
