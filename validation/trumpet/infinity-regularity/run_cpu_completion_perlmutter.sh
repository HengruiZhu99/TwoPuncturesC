#!/usr/bin/env bash
set -uo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
PUNCTURE_ROOT=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
bash "$PUNCTURE_ROOT/run_fd6_perlmutter.sh"
PUNCTURE_FD_STATUS=$?
bash "$PUNCTURE_ROOT/horizons_focus_perlmutter.sh"
PUNCTURE_HORIZON_STATUS=$?
printf '{"fd6_exit_code":%d,"horizon_exit_code":%d}\n' "$PUNCTURE_FD_STATUS" "$PUNCTURE_HORIZON_STATUS" > "$PUNCTURE_ROOT/focus-cpu-completion-status.json"
if [[ "$PUNCTURE_FD_STATUS" != 0 || "$PUNCTURE_HORIZON_STATUS" != 0 ]]; then exit 1;fi
