#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
cd "$base/validation-wave9999-v1"
export OMP_NUM_THREADS=16 OMP_PROC_BIND=spread OMP_PLACES=cores OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
python validation/wave22_final_check.py --library "$base/build-wave22-v1/libHiSpID.so" --checkpoint "$base/wave9999-study-v1/n160.hispid" --study "$base/wave9999-study-v1/result.json" --output "$base/wave9999-final-check-v1.json" > "$base/wave9999-final-check-v1.log" 2>&1
python validation/portable_sampler_migration.py --checkpoint "$base/wave9999-study-v1/n160.hispid" --producer-library "$base/build-wave22-v1/libHiSpID.so" --consumer-library "$base/build-wave22-cpu-v1/libHiSpID.so" --output "$base/wave9999-migration-v1.json" > "$base/wave9999-migration-v1.log" 2>&1
python "$base/athenak-wave22-v1/tst/test_suite/z4c/check_hispid_binary.py" --single --executable "$base/build-athenak-wave22-v1/src/athena" --checkpoint "$base/wave9999-study-v1/n160.hispid" --migration-proof "$base/wave9999-migration-v1.json" --output "$base/wave9999-horizon-v1" --levels 8,12,16 --geometry-threads 16 --strict-expansion --harmonic-storage factorized --guess-scale 1.00001 --domain-half-width 48 --reuse-shapes > "$base/wave9999-horizon-v1.log" 2>&1
