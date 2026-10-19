#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
cd "$base/validation-wave-radial-map-v2"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python validation/wave22_radial_profiles.py --library "$base/build-wave22-cpu-v1/libHiSpID.so" --checkpoints "$base/wave9999-study-v1/n96.hispid" "$base/wave9999-study-v1/n128.hispid" "$base/wave9999-study-v1/n160.hispid" --output "$base/wave9999-original-radial-profiles-v1" --points 257 --threads 8 --packet-constraints > "$base/wave9999-original-radial-profiles-v1.log" 2>&1
