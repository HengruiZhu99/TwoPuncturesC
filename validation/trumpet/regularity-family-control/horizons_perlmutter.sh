#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_DATA="$PUNCTURE_ROOT/qi-regularity-wide-v1/moderate224"
PUNCTURE_RUN="$PUNCTURE_ROOT/qi-regularity-wide-horizons-v1"
PUNCTURE_REUSE=${1:-fresh}
if [[ "$PUNCTURE_REUSE" == reuse-migration ]]; then
 PUNCTURE_RUN="$PUNCTURE_ROOT/qi-regularity-wide-horizons-v2"
elif [[ "$PUNCTURE_REUSE" != fresh ]]; then exit 2; fi
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
python3 - "$PUNCTURE_ROOT/qi-regularity-wide-v1/assessment.json" "$PUNCTURE_DATA/result.json" "$PUNCTURE_RUN/gate.json" <<'CHECK'
import json,sys
s=json.load(open(sys.argv[1]))
assert len(s['records'])==3 and s['increasing_resolutions']
assert all(x['physical_bounds_passed'] for x in s['records']),'all three absolute physical bounds must pass'
json.dump(dict(three_grid_passed=s['passed'],diagnostic_only=not s['passed'],
 reason='User requests finest AthenaK check; a failed residual-trend gate stays failed even if independent horizons pass'),open(sys.argv[3],'w'),indent=2)
r=json.load(open(sys.argv[2]));assert r['config']['seed_family']=='qi'
assert r['completed'] and r['diagnostics']['converged']
assert r['minimum_psi']>0 and r['minimum_metric_eigenvalue']>0 and r['exterior_stencil_unmodified']
for region in ('near','bulk'):
 q=r['physical'][region]
 assert q['H_rms']<1e-6 and q['M_rms']<1e-6 and q['H_max']<1e-4 and q['M_max']<1e-4
CHECK
printf '%s  %s\n' \
 c4c45eaee847a6ab17ad794c0d81321ea27d2a2cc313f84c8f0d622e7d8b4741 "$PUNCTURE_ROOT/build-arnoldi-20261018/libHiSpID.so" \
 14d0cf0ddc27fc64d1bed9fe6b4f8afbcce8f73ea73994af2ac17c28317dcc37 "$PUNCTURE_ROOT/build-axis-tau-polar-sampler/libHiSpID.so" \
 ff751123adcf4ea043c2b62d1f55ba10e5e966f462d7c606ef7b01e4decc57a5 "$PUNCTURE_ROOT/build-athenak-axis-tau-polar/src/athena" | sha256sum -c -
export PYTHONPATH="$PUNCTURE_ROOT/workflow-regularity-wide-20261018/python:$PUNCTURE_ROOT/workflow-regularity-wide-20261018/validation:$PUNCTURE_ROOT/workflow-regularity-wide-20261018/examples"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
if [[ "$PUNCTURE_REUSE" == reuse-migration ]]; then
 cp "$PUNCTURE_ROOT/qi-regularity-wide-horizons-v1/migration.json" "$PUNCTURE_RUN/migration.json"
else
 python3 "$PUNCTURE_ROOT/workflow-regularity-wide-20261018/validation/portable_sampler_migration.py" \
 --checkpoint "$PUNCTURE_DATA/diagnostic.checkpoint" \
 --producer-library "$PUNCTURE_ROOT/build-arnoldi-20261018/libHiSpID.so" \
 --consumer-library "$PUNCTURE_ROOT/build-axis-tau-polar-sampler/libHiSpID.so" \
 --output "$PUNCTURE_RUN/migration.json" > "$PUNCTURE_RUN/migration.log" 2>&1
fi
python3 "$PUNCTURE_ROOT/athenak-affine/tst/test_suite/z4c/check_hispid_binary.py" \
 --executable "$PUNCTURE_ROOT/build-athenak-axis-tau-polar/src/athena" \
 --checkpoint "$PUNCTURE_DATA/diagnostic.checkpoint" --allow-diagnostic \
 --migration-proof "$PUNCTURE_RUN/migration.json" --harmonic-storage factorized \
 --geometry-threads 16 --strict-expansion --reuse-shapes --levels 8,12,16 --timeout 1500 \
 --output "$PUNCTURE_RUN/surfaces" > "$PUNCTURE_RUN/horizons.log" 2>&1
