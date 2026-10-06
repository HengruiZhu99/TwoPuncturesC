#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_RUN="$PUNCTURE_ROOT/affine-horizon/exact"
mkdir -p "$PUNCTURE_RUN"
export PYTHONPATH="$PUNCTURE_ROOT/source-axisym/python:$PUNCTURE_ROOT/source-polar/examples:$PUNCTURE_ROOT/workflow-axisym/validation"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
python3 "$PUNCTURE_ROOT/source-polar/examples/export_athenak.py" --library "$PUNCTURE_ROOT/build-axisym-cuda/libHiSpID.so" --seed gamma10 --seed-family trumpet_r0_m --output "$PUNCTURE_RUN/gamma10.checkpoint" > "$PUNCTURE_RUN/export.json"
python3 "$PUNCTURE_ROOT/workflow-axisym/validation/portable_sampler_migration.py" --checkpoint "$PUNCTURE_RUN/gamma10.checkpoint" --producer-library "$PUNCTURE_ROOT/build-axisym-cuda/libHiSpID.so" --consumer-library "$PUNCTURE_ROOT/build-focused-sampler/libHiSpID.so" --output "$PUNCTURE_RUN/migration.json" > "$PUNCTURE_RUN/migration.log" 2>&1
python3 - "$PUNCTURE_RUN" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);d=json.loads((p/'export.json').read_text());d['migration_proof']=str(p/'migration.json');(p/'manifest.json').write_text(json.dumps(dict(gamma10=d),indent=2)+'\n')
PY
python3 "$PUNCTURE_ROOT/athenak-affine/tst/test_suite/z4c/check_hispid_controls.py" --executable "$PUNCTURE_ROOT/build-athenak-affine/src/athena" --manifest "$PUNCTURE_RUN/manifest.json" --cases gamma10 --affine-chart --boost-levels 4,8,12 --boost-flow-alpha .2 --flow-iterations 600 --worker-timeout 300 --harmonic-storage factorized --output "$PUNCTURE_RUN/horizons" > "$PUNCTURE_RUN/horizons.log" 2>&1
