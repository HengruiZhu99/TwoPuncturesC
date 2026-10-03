# Supplemental setup campaign

The separate setup matrix compares host and execution-space geometry in the
same new native image. Its126 workers use the declared40x80x16,80x160x16 and
128x256x28 grids, three repeats, Serial, OpenMP1/2/4/8/16 and one CUDA GPU.
It complements the frozen288-worker solve matrix and does not replace it.
No nonlinear or physical-binary acceptance follows from setup measurements.

Stage a committed validation source tree in a fresh scratch directory. Keep
the compiled Serial/OpenMP/CUDA images and their original source directories
unchanged. Write a JSON build map with exactly `serial`, `openmp`, and `cuda`
keys, whose values are those absolute build directories. Then run from the
staged validation tree with its `python`, `examples`, and `validation` on
`PYTHONPATH`:

```sh
python3 validation/kokkos_build_manifest.py \
  --root /absolute/original/scratch/root \
  --setup-build-map /absolute/setup-build-map.json \
  --output /absolute/setup-build-manifest.json
```

The builder reuses native/runtime dependency binding, records CMake caches,
checks the pinned Kokkos source, and compares the compiled native source
directory with the staged native source. Release, Kokkos, cubic row weighting
and benchmark instrumentation are required; explicit fast-math flags are
rejected. A manifest records compilation provenance and leaves actual setup
tests pending. Each variant must execute both strict setup and public seed
export suites before timings can qualify.

After the frozen solve campaign releases its allocation, request one GPU via
`shared_interactive`,32 logical CPUs, `gpu&hbm80g`, account `m3328_g`. Run the
benchmark sequentially through an `srun --gpus=1 --cpu-bind=cores` step:

```sh
python3 validation/benchmark_geometry_setup.py \
  --manifest /absolute/setup-build-manifest.json \
  --output validation/kokkos_setup_performance_20261003.json
```

For subsequent allocations use `--resume`, retaining the same manifest,
source, outputs and measured hardware class. Use the separately staged
`perlmutter_allocation_guard.py` with this results path to request a clean
between-worker checkpoint. The guard accepts the explicit126-worker inventory
as well as the legacy288-worker protocol; the frozen campaign retains its
original guard. `STOP_AFTER_WORKER` is honored before native controls and
timed workers. No allocation or numerical workers overlap.

Full source/artifact/control verification must succeed before completion or
qualified ratios are published. Failed suites, missing outputs, nonfinite
arrays or differences beyond the unchanged1e-10 preservation bound stay failed.
Retained timings may still characterize those attempted paths, with their
qualification state explicit.

RSS is a lifetime high-water through first sample, including input loading.
Kokkos callback peaks reset after initialization, retaining any tracked
allocations still live at reset. They exclude storage outside the callbacks,
including native vectors and CUDA stack/runtime. Driver observations are
whole-worker1Hz samples, including later verification, and give a lower bound
on transient VRAM peaks. Construction and first-sample timers are separate;
coefficient transformation is a subset of first sampling. The metadata-query
gap has its own wall-time field. Timing and memory definitions must remain
explicit in the final standalone LaTeX report.
