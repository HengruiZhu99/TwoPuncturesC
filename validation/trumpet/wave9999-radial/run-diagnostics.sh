#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
cd "$base/validation-wave-radial-map-v2"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
export OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
OMP_NUM_THREADS=4 python validation/wave22_radial_profiles.py --library "$base/build-wave-radial-map-cpu-v1/libHiSpID.so" --self-test --output "$base/wave9999-radial-map-oracle-v3.json" > "$base/wave9999-radial-map-oracle-v3.log" 2>&1
OMP_NUM_THREADS=16 python validation/wave22_residual_modes.py --library "$base/build-wave-radial-map-v1/libHiSpID.so" --checkpoint "$base/wave9999-radial-map-study-v6/n384.hispid" --output "$base/wave9999-radial-map-residual-modes-v1.json" > "$base/wave9999-radial-map-residual-modes-v1.log" 2>&1 &
mode_pid=$!
set +e
bash "$base/check-wave-original-radial-profiles-v1.sh"
profile_status=$?
wait "$mode_pid"
mode_status=$?
set -e
python - "$base" "$profile_status" "$mode_status" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);r=dict(legacy_profile_exit=int(sys.argv[2]),residual_localization_exit=int(sys.argv[3]),purpose='diagnostics, not acceptance of a refined solution')
(p/'wave9999-radial-diagnostics-exit.json').write_text(json.dumps(r,indent=2)+'\n')
print(r)
assert (p/'wave9999-original-radial-profiles-v1/result.json').exists(),'legacy diagnostic did not finish'
assert r['residual_localization_exit']==0,'residual localization failed'
PY
