#!/bin/bash
set -euo pipefail
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
salloc -A m3328_g -C gpu -q shared_interactive -N1 --gpus=1 --cpus-per-task=32 --mem=56G --ntasks=1 -t00:30:00 --job-name=hispid-half-turn srun --ntasks=1 --cpus-per-task=32 --gpus=1 bash "$base/run-half-turn-v2.sh"
