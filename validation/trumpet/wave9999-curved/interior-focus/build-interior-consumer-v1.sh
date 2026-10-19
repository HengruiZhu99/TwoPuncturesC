#!/bin/bash
set -euo pipefail
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
cmake -S "$base/source-wave-curved-v13" -B "$base/build-wave-curved-cpu-v13" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ -DPUNCTURES_KOKKOS=OFF -DHISPID_AXIS_TAU=ON -DHISPID_FOCUS_TAU=ON -DHISPID_FOCUS_ORDER=1 -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_STREAM_POLAR_PRECONDITIONER=ON -DHISPID_ROW_POWER=3 -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_BATCHED_MODAL_TRANSFER=ON -DHISPID_RADIAL_STRETCH=16 -DHISPID_ANGULAR_STRETCH=.3 -DHISPID_CURVED_MODAL_PRECONDITIONER=ON
cmake --build "$base/build-wave-curved-cpu-v13" --target HiSpID -j16
if [ ! -e "$base/athenak-curved-v13" ]; then
  cp -a "$base/athenak-wave22-v1" "$base/athenak-curved-v13"
  cp "$base/curved-hispid_checkpoint.hpp" "$base/athenak-curved-v13/src/pgen/z4c/hispid_checkpoint.hpp"
fi
bash "$base/build-athenak-curved-v13.sh"
