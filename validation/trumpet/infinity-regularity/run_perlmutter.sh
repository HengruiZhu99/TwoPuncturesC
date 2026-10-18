#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
module load cray-python/3.12.12
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
PUNCTURE_SOURCE=$PUNCTURE_ROOT/source-infinity2-20261018
PUNCTURE_BUILD=$PUNCTURE_ROOT/build-infinity2-20261018
PUNCTURE_RUN=$PUNCTURE_ROOT/infinity2-regularity-v2
mkdir "$PUNCTURE_RUN"
printf '%s\n' "$SLURM_JOB_ID" > "$PUNCTURE_RUN/job-id.txt"
cmake -S "$PUNCTURE_SOURCE" -B "$PUNCTURE_BUILD" -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper \
 -DPUNCTURES_KOKKOS=ON -DKokkos_DIR="$PUNCTURE_ROOT/build-cuda/cmake_packages/Kokkos" \
 -DHISPID_INFINITY_POWER=2 -DHISPID_AXIS_TAU=ON -DHISPID_EXACT_POLAR_TAU=ON \
 -DHISPID_ROW_POWER=3 -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_BATCHED_MODAL_TRANSFER=ON \
 > "$PUNCTURE_RUN/configure.log" 2>&1
cmake --build "$PUNCTURE_BUILD" -j 8 --target HiSpID test_hispid_axis test_hispid_kokkos > "$PUNCTURE_RUN/build.log" 2>&1
"$PUNCTURE_BUILD/test_hispid_axis" > "$PUNCTURE_RUN/axis-test.log" 2>&1
"$PUNCTURE_BUILD/test_hispid_kokkos" > "$PUNCTURE_RUN/kernel-test.log" 2>&1
sha256sum "$PUNCTURE_BUILD/"lib*.so > "$PUNCTURE_RUN/images.sha256"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores HISPID_TRACE_NEWTON=1
PUNCTURE_INITIAL=()
for PUNCTURE_N in 160 192 224; do
 python3 "$PUNCTURE_SOURCE/validation/run_trumpet_pilot.py" \
  --library "$PUNCTURE_BUILD/libHiSpID.so" --case qi-moderate-wide \
  --n "$PUNCTURE_N" --npolar "$((2*PUNCTURE_N))" --nphi 16 \
  --krylov-restart 64 --max-krylov 512 --linear-rtol .1 \
  --memory-mib 65536 --skip-observer-on-failure "${PUNCTURE_INITIAL[@]}" \
  --output "$PUNCTURE_RUN/moderate$PUNCTURE_N" > "$PUNCTURE_RUN/solve$PUNCTURE_N.log" 2>&1
 python3 - "$PUNCTURE_RUN/moderate$PUNCTURE_N/result.json" <<'CHECK'
import json,sys
r=json.load(open(sys.argv[1]));assert r['completed'] and r['diagnostics']['converged'],r['diagnostics']
assert r['minimum_psi']>0 and r['minimum_metric_eigenvalue']>0 and r['exterior_stencil_unmodified']
for region in ('near','bulk'):
 q=r['physical'][region];assert q['H_rms']<1e-6 and q['M_rms']<1e-6 and q['H_max']<1e-4 and q['M_max']<1e-4,(region,q)
CHECK
 PUNCTURE_INITIAL=(--initial "$PUNCTURE_RUN/moderate$PUNCTURE_N/diagnostic.checkpoint")
done
python3 "$PUNCTURE_SOURCE/validation/assess_trumpet_sequence.py" \
 --runs "$PUNCTURE_RUN/moderate160" "$PUNCTURE_RUN/moderate192" "$PUNCTURE_RUN/moderate224" \
 --output "$PUNCTURE_RUN/assessment.json" > "$PUNCTURE_RUN/assessment.log" 2>&1
