#!/usr/bin/env bash
# Sequential continuation of an already bound matrix. No builds, overlapping
# allocations, source/manifest edits, numerical retries or hardware-class waiver.
set -euo pipefail
export PUNCTURE_RUN_ROOT=${1:?isolated scratch run directory required}
PUNCTURE_MAX_ALLOCATIONS=${2:-8}
if [[ ! "$PUNCTURE_MAX_ALLOCATIONS" =~ ^[1-9][0-9]*$ ]]; then
  echo 'Positive maximum allocation count required' >&2; exit 2
fi
if [[ -n ${SLURM_JOB_ID:-} ]]; then
  echo 'Launch this observer from outside an allocation, after releasing the previous job' >&2; exit 2
fi
PUNCTURE_RESULTS="$PUNCTURE_RUN_ROOT/source/validation/kokkos_performance_20261002.json"
PUNCTURE_MANIFEST="$PUNCTURE_RUN_ROOT/performance-manifest.json"
PUNCTURE_STOP="$PUNCTURE_RUN_ROOT/source/validation/raw/kokkos_performance_20261002/STOP_AFTER_WORKER"
[[ -f "$PUNCTURE_RESULTS" && -f "$PUNCTURE_MANIFEST" && -f "$PUNCTURE_RUN_ROOT/allocation_guard.py" ]]
export PUNCTURE_RESULTS PUNCTURE_MANIFEST PUNCTURE_STOP
PUNCTURE_LOCK="$PUNCTURE_RUN_ROOT/performance_campaign.lock"
if ! mkdir -m700 "$PUNCTURE_LOCK"; then
  echo 'Campaign lock exists; inspect its owner and allocation before recovery' >&2; exit 2
fi
printf '%s %s\n' "$$" "$(hostname)" > "$PUNCTURE_LOCK/owner"
trap 'rm -f "$PUNCTURE_LOCK/owner"; rmdir "$PUNCTURE_LOCK"' EXIT
export PYTHONPATH="$PUNCTURE_RUN_ROOT/source/python:$PUNCTURE_RUN_ROOT/source/validation:$PUNCTURE_RUN_ROOT/source/examples"
export OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores
for ((PUNCTURE_CYCLE=1;PUNCTURE_CYCLE<=PUNCTURE_MAX_ALLOCATIONS;PUNCTURE_CYCLE++)); do
  if python3 -c 'import json,os,sys;sys.exit(0 if json.load(open(os.environ["PUNCTURE_RESULTS"])).get("declared_performance_completed") else 1)'; then
    echo 'Declared performance matrix completed; no further allocation requested'; exit 0
  fi
  PUNCTURE_PREVIOUS_JOB=$(python3 -c 'import json,os;print(json.load(open(os.environ["PUNCTURE_RESULTS"]))["allocation_epochs"][-1]["environment"]["SLURM_JOB_ID"])')
  PUNCTURE_PREVIOUS_JOBS=("$PUNCTURE_PREVIOUS_JOB")
  if [[ -f "$PUNCTURE_RUN_ROOT/performance_campaign_last_job.txt" ]]; then
    PUNCTURE_PREVIOUS_JOBS+=("$(cat "$PUNCTURE_RUN_ROOT/performance_campaign_last_job.txt")")
  fi
  for PUNCTURE_PREVIOUS_JOB in "${PUNCTURE_PREVIOUS_JOBS[@]}"; do
    [[ "$PUNCTURE_PREVIOUS_JOB" =~ ^[0-9]+$ ]]
    for ((PUNCTURE_WAIT=0;PUNCTURE_WAIT<120;PUNCTURE_WAIT++)); do
      PUNCTURE_PREVIOUS_STATE=$(squeue -h -j "$PUNCTURE_PREVIOUS_JOB" -o '%T')
      if [[ -z "$PUNCTURE_PREVIOUS_STATE" ]]; then break; fi
      if [[ "$PUNCTURE_PREVIOUS_STATE" != COMPLETING ]]; then
        echo "Previous allocation $PUNCTURE_PREVIOUS_JOB remains $PUNCTURE_PREVIOUS_STATE; no new allocation requested" >&2; exit 2
      fi
      sleep 1
    done
    if [[ -n "$PUNCTURE_PREVIOUS_STATE" ]]; then
      echo 'Previous allocation has not finished releasing resources' >&2; exit 2
    fi
  done
  echo "Requesting sequential shared allocation $PUNCTURE_CYCLE of $PUNCTURE_MAX_ALLOCATIONS"
  PUNCTURE_NODE_OPTIONS=()
  if [[ -n ${PUNCTURE_PREFERRED_NODE:-} ]]; then PUNCTURE_NODE_OPTIONS=(--nodelist="$PUNCTURE_PREFERRED_NODE"); fi
  salloc -N1 -q shared_interactive -C 'gpu&hbm80g' -G1 -n1 -c32 -t04:00:00 -A m3328_g -J codex-hispid-matrix --immediate=600 "${PUNCTURE_NODE_OPTIONS[@]}" bash --noprofile --norc -c '
    set -euo pipefail
    cd "$PUNCTURE_RUN_ROOT/source"
    printf "%s\n" "$SLURM_JOB_ID" > "$PUNCTURE_RUN_ROOT/performance_campaign_last_job.txt.tmp"
    mv "$PUNCTURE_RUN_ROOT/performance_campaign_last_job.txt.tmp" "$PUNCTURE_RUN_ROOT/performance_campaign_last_job.txt"
    rm -f "$PUNCTURE_STOP"
    python3 "$PUNCTURE_RUN_ROOT/allocation_guard.py" --job-id "$SLURM_JOB_ID" --results "$PUNCTURE_RESULTS" --wait-for-epoch > "$PUNCTURE_RUN_ROOT/guard-$SLURM_JOB_ID.log" 2>&1 &
    PUNCTURE_GUARD_PID=$!
    srun -n1 -c32 --gpus=1 --cpu-bind=cores python3 validation/benchmark_kokkos.py --manifest "$PUNCTURE_MANIFEST" --output "$PUNCTURE_RESULTS" --resume > "$PUNCTURE_RUN_ROOT/performance-$SLURM_JOB_ID.log" 2>&1 &
    PUNCTURE_WORKER_PID=$!
    trap "kill $PUNCTURE_GUARD_PID 2>/dev/null || true" EXIT
    while kill -0 "$PUNCTURE_WORKER_PID" 2>/dev/null && kill -0 "$PUNCTURE_GUARD_PID" 2>/dev/null; do sleep 1; done
    if ! kill -0 "$PUNCTURE_WORKER_PID" 2>/dev/null; then
      PUNCTURE_WORKER_STATUS=0
      wait "$PUNCTURE_WORKER_PID" || PUNCTURE_WORKER_STATUS=$?
      kill "$PUNCTURE_GUARD_PID" 2>/dev/null || true
      wait "$PUNCTURE_GUARD_PID" 2>/dev/null || true
      exit "$PUNCTURE_WORKER_STATUS"
    fi
    PUNCTURE_GUARD_STATUS=0
    wait "$PUNCTURE_GUARD_PID" || PUNCTURE_GUARD_STATUS=$?
    if ((PUNCTURE_GUARD_STATUS)); then touch "$PUNCTURE_STOP"; fi
    PUNCTURE_WORKER_STATUS=0
    wait "$PUNCTURE_WORKER_PID" || PUNCTURE_WORKER_STATUS=$?
    if ((PUNCTURE_GUARD_STATUS || PUNCTURE_WORKER_STATUS)); then
      echo "Observer/worker failed: $PUNCTURE_GUARD_STATUS/$PUNCTURE_WORKER_STATUS; retained logs and matrix are unchanged" >&2
      exit 1
    fi
  '
done
if python3 -c 'import json,os,sys;sys.exit(0 if json.load(open(os.environ["PUNCTURE_RESULTS"])).get("declared_performance_completed") else 1)'; then
  echo 'Declared performance matrix completed; no further allocation requested'; exit 0
fi
echo 'Allocation limit reached with matrix still incomplete; checkpoint retained' >&2
exit 3
