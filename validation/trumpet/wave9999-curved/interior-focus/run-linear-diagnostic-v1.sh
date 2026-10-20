#!/bin/bash
set -euo pipefail
module load python
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
cd "$base/source-wave-curved-v16"
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
export HISPID_TRACE_NEWTON=1 HISPID_PROBE_LINEAR_ACTION=1
python - "$base" <<'PY'
import json,sys,time
from pathlib import Path
from hispid import Backend
from checkpoint_export import read_checkpoint
base=Path(sys.argv[1]);c,u,meta=read_checkpoint(base/'wave9999-retained-block-seed-v1/n384.hispid')
c.max_newton=1;c.max_krylov=10;c.krylov_restart=24;c.memory_limit_mib=98304
b=Backend(str(base/'build-wave-curved-v16/libHiSpID.so'))
with b.create(c,execution='kokkos',wave=c.wave22,seed_family='liu') as s:
 s.set_unknowns(u);t=time.monotonic();d=s.solve(krylov='gmres',linear_rtol=.1)
 result=dict(purpose='bounded frozen physical linear diagnostic, not a converged solution',diagnostics=d,history=s.linear_history(),elapsed=time.monotonic()-t,source=meta,library_sha256=b.library_sha256())
 (base/'linear-diagnostic-v1.json').write_text(json.dumps(result,indent=2)+'\n')
 print(result,flush=True)
PY
