#!/usr/bin/env bash
set -euo pipefail
export PUNCTURE_RUN_ROOT=${1:?isolated scratch root required}
cd "$PUNCTURE_RUN_ROOT/source"
export PYTHONPATH=python:validation:examples OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores
python3 - <<'PY'
import json,os
from pathlib import Path
root=Path(os.environ['PUNCTURE_RUN_ROOT']);manifest=json.loads((root/'build-manifest.json').read_text())
manifest['variants']=[v for v in manifest['variants'] if v['id'] in ('reference','openmp16','cuda')]
(root/'preflight-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
PY
python3 validation/benchmark_kokkos.py --manifest "$PUNCTURE_RUN_ROOT/preflight-manifest.json" --grids 40:80:16 --repeats 1 --output validation/kokkos_preflight40_20261002.json
