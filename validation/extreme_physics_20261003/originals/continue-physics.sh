#!/usr/bin/env bash
set -euo pipefail
cd '/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/extreme-workflow-0606071'
export PYTHONPATH="$PWD/python:$PWD/examples:$PWD/validation"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=16 OMP_PROC_BIND=close OMP_PLACES=cores
PUNCTURE_EXTREME_ROOT='/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/extreme-human-override-20261003'
PUNCTURE_PHYSICS_LIBRARY='/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/build-cuda/libHiSpID.so'
PUNCTURE_PRESERVED_PERFORMANCE='/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/human-override-performance-preserved-20261003.json'
PUNCTURE_PERFORMANCE_ARTIFACTS='/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/source'
PUNCTURE_HUMAN_OVERRIDE="$PWD/validation/human_override_20261003.json"
PUNCTURE_PLAN="$PWD/validation/extreme_kokkos_plan.json"
run_floor() {
  local target=$1 separation=$2 grids=$3 status=0
  printf '%s\n' "source_floor_$target" > "$PUNCTURE_EXTREME_ROOT/phase.txt"
  srun -n1 -c32 --gpus=1 --cpu-bind=cores python3 validation/check_far_source_floor.py \
    --library "$PUNCTURE_PHYSICS_LIBRARY" --extreme --target-case "$target" \
    --execution kokkos --geometry host --threads 16 --memory-mib 32768 --seed-mass .5 \
    --coordinate-separation "$separation" --resolutions "$grids" \
    --output "$PUNCTURE_EXTREME_ROOT/source-floor-$target.json" \
    > "$PUNCTURE_EXTREME_ROOT/source-floor-$target.log" 2>&1 || status=$?
  printf '%s\n' "$status" > "$PUNCTURE_EXTREME_ROOT/source-floor-$target.exit"
  python3 -c 'import json,sys;d=json.load(open(sys.argv[1]));expected=2*len(sys.argv[2].split(","));sys.exit(0 if len(d.get("records",[]))==expected and all(r.get("completed") for r in d["records"]) else 1)' "$PUNCTURE_EXTREME_ROOT/source-floor-$target.json" "$grids"
}
run_case() {
  local case=$1 index=$2 target=$3 status=0
  printf '%s\n' "$case grid $index" > "$PUNCTURE_EXTREME_ROOT/phase.txt"
  srun -n1 -c32 --gpus=1 --cpu-bind=cores python3 validation/run_extreme_kokkos.py \
    --library "$PUNCTURE_PHYSICS_LIBRARY" --performance-results "$PUNCTURE_PRESERVED_PERFORMANCE" \
    --performance-artifact-root "$PUNCTURE_PERFORMANCE_ARTIFACTS" \
    --human-override "$PUNCTURE_HUMAN_OVERRIDE" --plan "$PUNCTURE_PLAN" \
    --seed-controls "$PUNCTURE_EXTREME_ROOT/seed-controls.json" \
    --source-floor-controls "$PUNCTURE_EXTREME_ROOT/source-floor-$target.json" \
    --case "$case" --grid-index "$index" --threads 16 --geometry host \
    --allow-diagnostic-investigation --output-directory "$PUNCTURE_EXTREME_ROOT/results" \
    > "$PUNCTURE_EXTREME_ROOT/$case-grid$index.log" 2>&1 || status=$?
  printf '%s\n' "$status" > "$PUNCTURE_EXTREME_ROOT/$case-grid$index.exit"
  python3 -c 'import json,sys;d=json.load(open(sys.argv[1]));index=int(sys.argv[3]);label=sys.argv[2]+("_fourier_control" if index==3 else "");rows=d.get(label,{}).get("records",[]);expected=1 if index==3 else index+1;sys.exit(0 if len(rows)>=expected else 1)' "$PUNCTURE_EXTREME_ROOT/results/results.json" "$case" "$index"
}
run_floor spin99 12 '80:160:16,104:208:20,128:256:28'
run_case aligned_spin99_kokkos 0 spin99
run_floor gamma10 25 '80:160:8,104:208:8,128:256:8,128:256:16'
run_case headon_gamma10_kokkos 0 gamma10
for index in 1 2; do
  run_case aligned_spin99_kokkos "$index" spin99
  run_case headon_gamma10_kokkos "$index" gamma10
done
run_case headon_gamma10_kokkos 3 gamma10
printf '%s\n' 'initial_binary_grids_retained' > "$PUNCTURE_EXTREME_ROOT/phase.txt"
