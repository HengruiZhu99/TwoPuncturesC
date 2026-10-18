#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-final-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/final-order-control-v1
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
cmake -S "$PUNCTURE_SOURCE" -B "$PUNCTURE_ROOT/build-final-c1-sampler" \
 -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ -DCMAKE_BUILD_TYPE=Release \
 -DHISPID_AXIS_TAU=ON -DHISPID_FOCUS_TAU=ON -DHISPID_FOCUS_ORDER=1 -DHISPID_EXACT_POLAR_TAU=ON \
 -DHISPID_ROW_POWER=3 -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_BATCHED_MODAL_TRANSFER=ON > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_ROOT/build-final-c1-sampler" -j8 --target HiSpID test_hispid_focus > "$PUNCTURE_RUN/build.log" 2>&1
"$PUNCTURE_ROOT/build-final-c1-sampler/test_hispid_focus" > "$PUNCTURE_RUN/focus-unit.log"
/opt/cray/pe/gcc-native/14/bin/g++ -O2 -std=c++17 -I"$PUNCTURE_SOURCE/src" -I"$PUNCTURE_SOURCE/include" "$PUNCTURE_SOURCE/validation/focus-order/check_focus_order_compatibility.cpp" -o "$PUNCTURE_RUN/weights-check"
python3 "$PUNCTURE_SOURCE/validation/focus-order/check_focus_order_compatibility.py" --executable "$PUNCTURE_RUN/weights-check" --runs "$PUNCTURE_ROOT/focus-tau-v1" --output "$PUNCTURE_RUN/compatibility" > "$PUNCTURE_RUN/compatibility.log" 2>&1
sha256sum "$PUNCTURE_ROOT/build-final-c1-sampler/"lib*.so "$PUNCTURE_RUN/weights-check" > "$PUNCTURE_RUN/images.sha256"
for PUNCTURE_WAIT in $(seq 1 40); do
 if python3 - "$PUNCTURE_ROOT/focus2-refine-v1/moderate192/result.json" <<'READY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);raise SystemExit(0 if p.exists() and json.loads(p.read_text()).get('completed') else 1)
READY
 then break;fi
 sleep 10
done
python3 - "$PUNCTURE_ROOT/focus2-refine-v1/moderate192/result.json" "$PUNCTURE_RUN/selection.json" <<'SELECT'
import json,sys
r=json.load(open(sys.argv[1]));assert r['completed']
passed=bool(r['diagnostics']['converged'] and all(q['H_rms']<1e-6 and q['M_rms']<1e-6 and q['H_max']<1e-4 and q['M_max']<1e-4 for k,q in r['physical'].items() if k in ('near','bulk')) and r['physical']['bulk']['M_rms']<1.0538120185710488e-6/4 and r['physical']['bulk']['H_rms']<1.9591214219981626e-9)
json.dump(dict(c2_prefix_passed=passed,selected_for_horizon='c2_pending_sequence' if passed else 'c1_finest224',reason='Predeclared C2 refinement gate; original physical acceptance remains separately assessed.'),open(sys.argv[2],'w'),indent=2)
SELECT
if python3 - "$PUNCTURE_RUN/selection.json" <<'SELECT'
import json,sys
raise SystemExit(1 if json.load(open(sys.argv[1]))['c2_prefix_passed'] else 0)
SELECT
then
 bash "$PUNCTURE_ROOT/horizons_focus_c1_final.sh"
fi
