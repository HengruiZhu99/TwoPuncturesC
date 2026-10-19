#!/bin/bash
set -euo pipefail
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
# Wait for the current bounded pilot's completed report, written AFTER its
# checkpoint. It is an initial guess even if its native gate still fails.
while true; do
 if python3 - "$base/wave9999-retained-block-seed-v1/result.json" <<'PY'
from pathlib import Path
import json,sys
try:
 p=Path(sys.argv[1]);d=json.loads(p.read_text());r=d['levels'][0]
 assert len(d['levels'])==1 and r['n']==[384,128,16] and 'checkpoint' in r
 assert (p.parent/'n384.hispid').is_file()
 print('Saved initial-guess checkpoint; native status:',r['diagnostics'],flush=True)
except (OSError,ValueError,KeyError,IndexError,AssertionError):raise SystemExit(1)
PY
 then break; fi
 if ! squeue -h -j 60008107 2>/dev/null | grep -q .; then echo 'Pilot ended without a completed checkpoint receipt'; exit 1; fi
 sleep 10
done
# No overlapping GPU allocations. Preserve the completed pilot before the
# smaller-restart full study can consume further compute time.
scancel 60008107 || true
while squeue -h -j 60008107 2>/dev/null | grep -q .; do sleep 2; done
salloc -A m3328_g -C gpu -q shared_interactive -N1 --gpus=1 --cpus-per-task=32 --ntasks=1 -t02:00:00 --job-name=hispid-outer24 srun --ntasks=1 --cpus-per-task=32 --gpus=1 bash "$base/run-large-outer-all-v1.sh"
