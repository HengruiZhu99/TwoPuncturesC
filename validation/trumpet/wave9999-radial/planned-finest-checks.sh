#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
cd "$base/validation-wave-radial-map-v2"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python validation/wave22_radial_profiles.py --library "$base/build-wave-radial-map-cpu-v1/libHiSpID.so" --checkpoints "$base/wave9999-radial-map-study-v6/n384.hispid" "$base/wave9999-radial-map-study-v6/n480.hispid" "$base/wave9999-radial-map-study-v6/n576.hispid" --output "$base/wave9999-radial-map-profiles-v6" --points 257 --threads 16 --packet-constraints > "$base/wave9999-radial-map-profiles-v6.log" 2>&1
python validation/portable_sampler_migration.py --checkpoint "$base/wave9999-radial-map-study-v6/n576.hispid" --producer-library "$base/build-wave-radial-map-v1/libHiSpID.so" --consumer-library "$base/build-wave-radial-map-cpu-v1/libHiSpID.so" --output "$base/wave9999-radial-map-migration-v6.json" > "$base/wave9999-radial-map-migration-v6.log" 2>&1
python - <<'PY'
from pathlib import Path
import numpy as np
base=Path('/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005')
coeff=np.loadtxt(base/'wave9999-horizon-v1/quadrature_l16_n48/hispid.horizon_shape_0.txt')
np.savetxt(base/'wave9999-radial-map-guess-v6.txt',coeff[:81][None,:],fmt='%.17g')
PY
LD_LIBRARY_PATH="$base/build-wave-radial-map-cpu-v1:${LD_LIBRARY_PATH:-}" python "$base/athenak-wave22-v1/tst/test_suite/z4c/check_hispid_binary.py" --single --executable "$base/build-athenak-wave22-v1/src/athena" --checkpoint "$base/wave9999-radial-map-study-v6/n576.hispid" --migration-proof "$base/wave9999-radial-map-migration-v6.json" --output "$base/wave9999-radial-map-horizon-v6" --levels 8,12,16 --geometry-threads 16 --strict-expansion --harmonic-storage factorized --initial-shapes "$base/wave9999-radial-map-guess-v6.txt" --domain-half-width 48 --reuse-shapes > "$base/wave9999-radial-map-horizon-v6.log" 2>&1
