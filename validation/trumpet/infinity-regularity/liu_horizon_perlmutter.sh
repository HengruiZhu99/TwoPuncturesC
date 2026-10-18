#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-infinity2-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/liu-exact-horizon-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export PYTHONPATH="$PUNCTURE_SOURCE/python:$PUNCTURE_SOURCE/validation:$PUNCTURE_SOURCE/examples"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
python3 "$PUNCTURE_SOURCE/examples/export_athenak.py" \
 --library "$PUNCTURE_ROOT/build-infinity2-20261018/libHiSpID.so" \
 --seed kerr99 --seed-family liu --output "$PUNCTURE_RUN/kerr99.checkpoint" > "$PUNCTURE_RUN/export.json"
python3 "$PUNCTURE_SOURCE/validation/portable_sampler_migration.py" \
 --checkpoint "$PUNCTURE_RUN/kerr99.checkpoint" \
 --producer-library "$PUNCTURE_ROOT/build-infinity2-20261018/libHiSpID.so" \
 --consumer-library "$PUNCTURE_ROOT/build-infinity2-sampler/libHiSpID.so" \
 --output "$PUNCTURE_RUN/migration.json" > "$PUNCTURE_RUN/migration.log" 2>&1
python3 - "$PUNCTURE_RUN" <<'PY'
import json,sys
from pathlib import Path
r=Path(sys.argv[1]);entry=json.loads((r/'export.json').read_text());entry['migration_proof']=str(r/'migration.json')
(r/'manifest.json').write_text(json.dumps(dict(kerr99=entry),indent=2)+'\n')
PY
python3 "$PUNCTURE_ROOT/athenak-infinity2-20261018/tst/test_suite/z4c/check_hispid_controls.py" \
 --executable "$PUNCTURE_ROOT/build-athenak-infinity2/src/athena" \
 --manifest "$PUNCTURE_RUN/manifest.json" --cases kerr99 --harmonic-storage factorized \
 --output "$PUNCTURE_RUN/surfaces" > "$PUNCTURE_RUN/horizons.log" 2>&1
