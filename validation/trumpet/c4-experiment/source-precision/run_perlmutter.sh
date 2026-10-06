#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?compute allocation required}"
root=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
out="$root/c4-source-precision-v1"
mkdir "$out"
printf '%s\n' "$SLURM_JOB_ID" > "$out/job-id.txt"
/opt/cray/pe/gcc-native/14/bin/g++ -std=c++17 -O2 -I"$root/source-c4-backtrack/src" -I"$root/source-c4-backtrack/include" "$root/diagnose_trumpet_source_precision.cpp" -o "$out/probe"
sha256sum "$root/diagnose_trumpet_source_precision.cpp" "$root/source-c4-backtrack/src/HiSpID_geometry_kernels.hpp" "$root/source-c4-backtrack/src/HiSpID_cache_kernels.hpp" "$out/probe" > "$out/hashes.txt"
"$out/probe" > "$out/source.csv"
