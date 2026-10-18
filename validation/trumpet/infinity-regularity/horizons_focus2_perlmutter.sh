#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_DATA="$PUNCTURE_ROOT/focus2-v2/moderate224"
PUNCTURE_RUN="$PUNCTURE_ROOT/focus2-horizons-v1"
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
python3 - "$PUNCTURE_ROOT/focus2-v2/assessment.json" "$PUNCTURE_DATA/result.json" "$PUNCTURE_RUN/gate.json" <<'CHECK'
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
sha256sum -c "$PUNCTURE_ROOT/focus2-consumer-v1/images.sha256"
sha256sum -c "$PUNCTURE_ROOT/focus2-build-v3/images.sha256"
export PYTHONPATH="$PUNCTURE_ROOT/source-focus2-20261018/python:$PUNCTURE_ROOT/source-focus2-20261018/validation:$PUNCTURE_ROOT/source-focus2-20261018/examples"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python3 "$PUNCTURE_ROOT/source-focus2-20261018/validation/portable_sampler_migration.py" \
 --checkpoint "$PUNCTURE_DATA/diagnostic.checkpoint" \
 --producer-library "$PUNCTURE_ROOT/build-focus2-cuda/libHiSpID.so" \
 --consumer-library "$PUNCTURE_ROOT/build-focus2-sampler/libHiSpID.so" \
 --output "$PUNCTURE_RUN/migration.json" > "$PUNCTURE_RUN/migration.log" 2>&1
python3 "$PUNCTURE_ROOT/athenak-focus2-20261018/tst/test_suite/z4c/check_hispid_binary.py" \
 --executable "$PUNCTURE_ROOT/build-athenak-focus2/src/athena" \
 --checkpoint "$PUNCTURE_DATA/diagnostic.checkpoint" --allow-diagnostic \
 --migration-proof "$PUNCTURE_RUN/migration.json" --harmonic-storage factorized \
 --geometry-threads 16 --strict-expansion --reuse-shapes --levels 8,12,16 --timeout 1500 \
 --output "$PUNCTURE_RUN/surfaces" > "$PUNCTURE_RUN/horizons.log" 2>&1
