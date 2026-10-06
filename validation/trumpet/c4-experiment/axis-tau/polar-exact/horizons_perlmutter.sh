#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_DATA="$PUNCTURE_ROOT/axis-tau-polar-sequence-v1/moderate256"
PUNCTURE_RUN="$PUNCTURE_ROOT/axis-tau-polar-horizon-v1"
# Individual physical bounds gate this consumer run; the separate three-grid
# assessment still governs convergence and is never inherited by this test.
python3 - "$PUNCTURE_DATA/result.json" <<'CHECK'
import json,sys
x=json.load(open(sys.argv[1]))
assert x['completed'] and x['diagnostics']['converged']
assert x['minimum_psi']>0 and x['minimum_metric_eigenvalue']>0 and x['exterior_stencil_unmodified']
for region in ('near','bulk'):
 p=x['physical'][region]
 assert p['H_rms']<1e-6 and p['M_rms']<1e-6 and p['H_max']<1e-4 and p['M_max']<1e-4
CHECK
printf '%s  %s\n' \
 9489bca44c0bbd6ea20ae4996b4e2cdf225e84e3436ec80b253c4ffaf2b0c3e6 "$PUNCTURE_ROOT/build-axis-tau-polar/libHiSpID.so" \
 14d0cf0ddc27fc64d1bed9fe6b4f8afbcce8f73ea73994af2ac17c28317dcc37 "$PUNCTURE_ROOT/build-axis-tau-polar-sampler/libHiSpID.so" \
 ff751123adcf4ea043c2b62d1f55ba10e5e966f462d7c606ef7b01e4decc57a5 "$PUNCTURE_ROOT/build-athenak-axis-tau-polar/src/athena" | sha256sum -c -
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export PYTHONPATH="$PUNCTURE_ROOT/source-axis-tau-polar/python:$PUNCTURE_ROOT/source-axis-tau-polar/validation:$PUNCTURE_ROOT/source-axis-tau-polar/examples"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python3 "$PUNCTURE_ROOT/source-axis-tau-polar/validation/portable_sampler_migration.py" \
 --checkpoint "$PUNCTURE_DATA/diagnostic.checkpoint" \
 --producer-library "$PUNCTURE_ROOT/build-axis-tau-polar/libHiSpID.so" \
 --consumer-library "$PUNCTURE_ROOT/build-axis-tau-polar-sampler/libHiSpID.so" \
 --output "$PUNCTURE_RUN/migration.json" > "$PUNCTURE_RUN/migration.log" 2>&1
python3 "$PUNCTURE_ROOT/athenak-affine/tst/test_suite/z4c/check_hispid_binary.py" \
 --executable "$PUNCTURE_ROOT/build-athenak-axis-tau-polar/src/athena" \
 --checkpoint "$PUNCTURE_DATA/diagnostic.checkpoint" --allow-diagnostic \
 --migration-proof "$PUNCTURE_RUN/migration.json" --harmonic-storage factorized \
 --geometry-threads 16 --strict-expansion --levels 8,12,16 --timeout 1500 \
 --output "$PUNCTURE_RUN/surfaces" > "$PUNCTURE_RUN/horizons.log" 2>&1
