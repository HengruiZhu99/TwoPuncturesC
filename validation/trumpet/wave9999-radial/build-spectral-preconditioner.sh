#!/bin/bash
set -euo pipefail
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export NVCC_WRAPPER_DEFAULT_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++
cmake -S "$base/source-wave-radial-map-v1" -B "$base/build-wave-radial-spectral-v1" -DCMAKE_BUILD_TYPE=Release -DHISPID_AXIS_TAU=ON -DHISPID_FOCUS_TAU=ON -DHISPID_FOCUS_ORDER=1 -DHISPID_EXACT_POLAR_TAU=ON -DHISPID_ROW_POWER=3 -DHISPID_COMPENSATED_MODAL_PROJECTION=ON -DHISPID_BATCHED_MODAL_TRANSFER=ON -DHISPID_RADIAL_STRETCH=1 -DHISPID_SPECTRAL_RADIAL_PRECONDITIONER=ON -DPUNCTURES_KOKKOS=ON -DCMAKE_CXX_COMPILER=/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/kokkos/bin/nvcc_wrapper -DKokkos_DIR="$base/build-cuda/cmake_packages/Kokkos"
cmake --build "$base/build-wave-radial-spectral-v1" --target HiSpID test_hispid_radial_preconditioner -j16
"$base/build-wave-radial-spectral-v1/test_hispid_radial_preconditioner"
sha256sum "$base/build-wave-radial-spectral-v1/libHiSpID.so"
