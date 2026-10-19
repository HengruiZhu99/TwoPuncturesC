#!/bin/bash
set -euo pipefail
base=/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
cmake -S "$base/athenak-curved-v13" -B "$base/build-athenak-curved-v13" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/opt/cray/pe/gcc-native/14/bin/g++ -DPROBLEM=z4c/hispid -DAthena_ENABLE_OPENMP=ON -DAthena_ENABLE_MPI=OFF -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_SERIAL=OFF -DKokkos_ENABLE_CUDA=OFF -DHISPID_ROOT="$base/source-wave-curved-v13" -DHISPID_LIBRARY_DIR="$base/build-wave-curved-cpu-v13"
cmake --build "$base/build-athenak-curved-v13" -j16
sha256sum "$base/build-athenak-curved-v13/src/athena" "$base/build-wave-curved-cpu-v13/libHiSpID.so"
