#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
study=${1:?accepted three-grid study directory required}
producer=${2:?exact producer library required}
tag=${3:?new unique artifact tag required}
cd "$base/post-wave-curved-v13"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python validation/wave22_near_hole.py --library "$base/build-wave-curved-cpu-v13/libHiSpID.so" --checkpoints "$study/n384.hispid" "$study/n480.hispid" "$study/n576.hispid" --output "$base/$tag-near-hole" --threads 16
python validation/wave22_radial_profiles.py --library "$base/build-wave-curved-cpu-v13/libHiSpID.so" --checkpoints "$study/n384.hispid" "$study/n480.hispid" "$study/n576.hispid" --output "$base/$tag-profiles" --points 257 --threads 16 --packet-constraints
python validation/portable_sampler_migration.py --checkpoint "$study/n576.hispid" --producer-library "$producer" --consumer-library "$base/build-wave-curved-cpu-v13/libHiSpID.so" --output "$base/$tag-migration.json"
python - "$base/$tag-guess.txt" <<'PY'
from pathlib import Path
import numpy as np,sys
base=Path('/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005')
coeff=np.loadtxt(base/'wave9999-horizon-v1/quadrature_l16_n48/hispid.horizon_shape_0.txt')
np.savetxt(sys.argv[1],coeff[:81][None,:],fmt='%.17g')
PY
LD_LIBRARY_PATH="$base/build-wave-curved-cpu-v13:${LD_LIBRARY_PATH:-}" python "$base/athenak-curved-v13/tst/test_suite/z4c/check_hispid_binary.py" --single --executable "$base/build-athenak-curved-v13/src/athena" --checkpoint "$study/n576.hispid" --migration-proof "$base/$tag-migration.json" --output "$base/$tag-horizon" --levels 8,12,16 --geometry-threads 16 --strict-expansion --harmonic-storage factorized --initial-shapes "$base/$tag-guess.txt" --domain-half-width 48 --reuse-shapes
