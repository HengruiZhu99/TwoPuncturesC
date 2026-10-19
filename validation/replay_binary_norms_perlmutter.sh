#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-low-single-final-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/lower-resolution-20261018-v1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=16
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/norm-job-id.txt"
for PUNCTURE_CASE in moderate moderate-prefix paper-binary; do
 python3 "$PUNCTURE_SOURCE/validation/binary_reference_norm.py" \
  --study "$PUNCTURE_RUN/$PUNCTURE_CASE" \
  --by-library "$PUNCTURE_ROOT/build-by-capacity-20261018/libTwoPunctures.so" \
  --hispid-library "$PUNCTURE_ROOT/build-focus-dd-cuda/libHiSpID.so" \
  --output "$PUNCTURE_RUN/$PUNCTURE_CASE-reference-norm-bound.json"
done
python3 - "$PUNCTURE_RUN" <<'PY'
import json,os,platform,sys,numpy,scipy
from pathlib import Path
p=Path(sys.argv[1]); (p/'environment.json').write_text(json.dumps(dict(
 solve_job='59970516', saved_norm_job=os.environ['SLURM_JOB_ID'],
 python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__, platform=platform.platform(),
 scope='one shared interactive A100,32 allocated CPU cores; sequential native GPU solves,16 host threads; direct spherical solves and saved-norm replay on CPU'),indent=2)+'\n')
PY
