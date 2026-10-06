#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?allocated compute node required}"
PUNCTURE_ROOT=${1:?run root}
export PYTHONPATH="$PUNCTURE_ROOT/source/python:$PUNCTURE_ROOT/source/examples:$PUNCTURE_ROOT/source/validation"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=false
PUNCTURE_OUT="$PUNCTURE_ROOT/first-physics"
mkdir "$PUNCTURE_OUT"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_OUT/job-id.txt"
python3 "$PUNCTURE_ROOT/source/examples/export_athenak.py" --library "$PUNCTURE_ROOT/build-cuda/libHiSpID.so" --seed kerr99 --seed-family trumpet_r0_m --output "$PUNCTURE_OUT/kerr99.checkpoint" > "$PUNCTURE_OUT/kerr99-export.json"
python3 "$PUNCTURE_ROOT/source/validation/portable_sampler_migration.py" --checkpoint "$PUNCTURE_OUT/kerr99.checkpoint" --producer-library "$PUNCTURE_ROOT/build-cuda/libHiSpID.so" --consumer-library "$PUNCTURE_ROOT/build-sampler/libHiSpID.so" --output "$PUNCTURE_OUT/kerr99-migration.json" > "$PUNCTURE_OUT/migration.log" 2>&1
python3 - "$PUNCTURE_OUT" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);entry=json.loads((p/'kerr99-export.json').read_text());entry['migration_proof']=str(p/'kerr99-migration.json');(p/'manifest.json').write_text(json.dumps(dict(kerr99=entry),indent=2)+'\n')
PY
# A failed seed search is retained; the independent moderate pilot may diagnose
# elliptic behavior, but neither activity grants physical binary acceptance.
set +e
python3 "$PUNCTURE_ROOT/athenak/tst/test_suite/z4c/check_hispid_controls.py" --executable "$PUNCTURE_ROOT/build-athenak/src/athena" --manifest "$PUNCTURE_OUT/manifest.json" --output "$PUNCTURE_OUT/horizons" --cases kerr99 --harmonic-storage factorized > "$PUNCTURE_OUT/horizons.log" 2>&1
PUNCTURE_HORIZON_STATUS=$?
printf '%s\n' "$PUNCTURE_HORIZON_STATUS" > "$PUNCTURE_OUT/horizon-exit.txt"
set -e
python3 "$PUNCTURE_ROOT/source/validation/run_trumpet_pilot.py" --library "$PUNCTURE_ROOT/build-cuda/libHiSpID.so" --output "$PUNCTURE_OUT/moderate-pilot" > "$PUNCTURE_OUT/pilot.log" 2>&1
printf 'completed\n' > "$PUNCTURE_OUT/completed.txt"
