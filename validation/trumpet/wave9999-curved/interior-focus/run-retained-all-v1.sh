#!/bin/bash
set -euo pipefail
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
bash "$base/run-retained-block-v1.sh"
bash "$base/run-retained-post-v1.sh" "$base/wave9999-retained-block-study-v1" "$base/build-wave-curved-v16/libHiSpID.so" wave9999-retained-block-final
