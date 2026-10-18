#!/usr/bin/env bash
set -euo pipefail
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
printf '{"deferred":true,"reason":"Both fixed three-grid observers fail bulk trend. Preserve this prototype; measure horizons on the final convergence candidate."}\n' > "$PUNCTURE_ROOT/focus-tau-horizon-deferred.json"
python3 "$PUNCTURE_ROOT/source-higher-observer-20261018/validation/diagnose_focus_c2.py" --runs "$PUNCTURE_ROOT/focus-tau-v1" --output "$PUNCTURE_ROOT/focus-c2-compatibility.json"
