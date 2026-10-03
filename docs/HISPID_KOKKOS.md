# Explicit Kokkos execution

The optional CMake build provides Serial, OpenMP or CUDA execution for either
Bowen–York or HiSpID and either GMRES or BiCGStab. The execution space is fixed
by the image; the equation system, Krylov method and execution choice are
separate options. The historical Makefile and reference entry points remain
the default. Actual execution-space kernel controls, original BY fixtures and
bitwise reference-default controls pass. HiSpID's strict raw-P/coefficient gate
remains failed. The full performance matrix is running; its report will
precede the new physical studies.

Both execution policies instantiate one C-compatible Krylov controller.
Reference vector operations preserve the original C arithmetic and BY's
historical reductions. The Kokkos policy keeps all inner-solve vectors,
reductions, spectral derivatives and modal preconditioner applications in the
chosen memory space. An inherited host reduction callback is explicitly
replaced in the BY adapter; arbitrary host callbacks are rejected by the
generic device solver. True-residual recomputation, restarts, limits and
failure handling use the same controller.

HiSpID solves four coupled equations: one Hamiltonian and three momentum
constraints. Its seed conformal metric, mean curvature and lapse prescription
are free data, not additional elliptic unknowns. Geometry uses the existing
high-precision CPU construction. Newton control, line search, nonlinear host
iterates, native block factorization and physical sampling remain host work.
The OpenMP build parallelizes geometry and independent modal factor groups;
CUDA transfers fixed geometry once and keeps the inner iteration resident.
The reference BY nonlinear residual and Newton loop retain their C arithmetic.
The opt-in BY path uses one scoped Kokkos workspace throughout each Newton
call, reusing differentiation matrices and fixed seed data for nonlinear
residuals and linear JVPs. It exports every raw AB/phi and physical Cartesian
derivative array expected by the native solver. Mass changes refresh the
potential; changed separation, momenta or spins invalidate the workspace.
The residual-debug-file option is rejected for this explicit execution path.
Native final-state sampling and modal factorization remain CPU work.
Kokkos residual convergence is a candidate stopping point: the original
`F_of_v` confirms the same outer tolerance. If confirmation fails, remaining
outer residuals use the original routine for polishing with the resident
Kokkos linear action. The tolerance, forcing and Newton limit are unchanged;
confirmation calls, polishing steps and both norms are recorded separately.

The modal adapter shares factors between Fourier partners and, for HiSpID,
between the three vector components. Host execution borrows native factors;
it does not pack duplicate LU/transfer banks or retain replicated FD stencils.
Unused reference derivative banks are omitted in HiSpID Kokkos contexts, and
BY does not allocate the regular-mode scratch used only by HiSpID.
CUDA imports factors and inverts each pivoted radial Schur block once per
preconditioner setup for row-parallel products. Existing dense transfer/lower
banks are transposed in place for coalesced CUDA row access, without duplicate
storage. This opt-in arithmetic requires
its own matrix-action and complete-state qualification. It does not alter the
reference triangular solves.

## Build

Use GSL and a pinned Kokkos4.7.2 source checkout (commit
6739bc623081648af9e752b616d9671527922cbf), or a shared Kokkos package.
Keep each build and native image in a distinct directory. A shared Kokkos
runtime is required because HiSpID and TwoPunctures are separate DSOs.

```sh
cmake -S . -B build-kokkos-openmp \
  -DCMAKE_BUILD_TYPE=Release -DPUNCTURES_KOKKOS=ON \
  -DPUNCTURES_KOKKOS_SOURCE=/absolute/path/to/kokkos \
  -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_SERIAL=ON
cmake --build build-kokkos-openmp -j1
OMP_NUM_THREADS=1 ctest --test-dir build-kokkos-openmp --output-on-failure -j1
```

For one Perlmutter A100, use nvcc_wrapper as C++ compiler and add
`Kokkos_ENABLE_CUDA=ON`, `Kokkos_ENABLE_CUDA_LAMBDA=ON`,
`Kokkos_ARCH_AMPERE80=ON`, and host OpenMP. Native C remains C99/GCC.
The performance script uses cubic row weighting (`HISPID_ROW_POWER=3`)
and timing-only BY instrumentation (`PUNCTURES_BENCHMARK=ON`); ordinary
builds retain the default sixth-power HiSpID weighting.

Request `shared_interactive`, one GPU,32 logical CPUs (16 physical cores) and
account m3328_g, then run every build/test/benchmark through `srun --gpus=1`.
The isolated `validation/perlmutter_kokkos_build.sh` builds all images
sequentially. Do not reuse a full-node allocation for a one-GPU study.
For a resumed performance allocation use `-C "gpu&hbm80g"` to retain the
measured 80 GB GPU class; the coordinator also verifies CPU/cache topology,
affinity count and GPU memory/driver class. Keep its frozen source/images and
manifest unchanged throughout the matrix.

An administrative observer can request a clean checkpoint before Slurm's
deadline. Copy `validation/perlmutter_allocation_guard.py` outside the frozen
source tree and run it with the current `--job-id` and absolute `--results`
path. It sets `STOP_AFTER_WORKER` with enough time for the full worker timeout
and finalization margin, without killing work. Release that allocation after
the coordinator exits; request another one-GPU shared allocation, remove the
marker, and restart the same command with `--resume`. The observer exits if
the allocation epoch changes. No numerical workers may overlap.

## Python and ownership

```python
from hispid import Backend
from execution import select
backend = Backend('/absolute/build-kokkos-openmp/libHiSpID.so')
select(backend.lib, 'kokkos', threads=4)
config = backend.config()
with backend.create(config, execution='kokkos') as data:
    diagnostics = data.solve(linear_rtol=.001, krylov='gmres')
```

BY's `Solution.solve` accepts `execution='kokkos'` with
`preconditioner='modal'` and an explicit positive `linear_rtol`. It supports
both Krylov methods. Unavailable execution spaces fail explicitly.
Distinct builds are compared in fresh processes; adapters verify the primary
and dependent puncture symbol images, and record the mapped Kokkos runtime.

Kokkos initialization is serialized and occurs once. Explicit repeated thread
requests must match actual host concurrency. An externally owned runtime is
never finalized or given replacement allocation callbacks. Owned initialization
registers teardown after Kokkos's static state. Numerical Kokkos operations,
initialization and statistics access are serialized process-wide. The caller
must not mutate or destroy a context while another caller uses it. Native BY's
global parameter setup remains non-reentrant and its Python adapter holds a
process-wide context lock.

HiSpID's explicit memory budget maximum is extended8192 to65536MiB; the default
2048MiB and the existing reference estimate remain unchanged. Kokkos adds a
conservative aggregate bound, including maximum lazy Krylov storage and
factor/inverse overlap, before context allocation. BY uses
`TP_execution_memory_limit_mib` (default8192) for its aggregate estimate.
CUDA also checks available device memory with512MiB headroom; a failed CUDA
memory query rejects allocation. Estimates are conservative capacity screens,
not memory measurements.

Allocation callbacks report Kokkos host/device allocations observed after
runtime initialization, including current/peak bytes and setup overlap. They exclude native std::vector/GSL storage and CUDA runtime
allocations. Pinned memory is host; unified/managed memory is classified by its
logical Kokkos device space. Host RSS measures the whole worker. Driver process
memory is sampled separately at1Hz and cannot resolve every transient peak.
Logical resident/workspace estimates and API copy-byte counters are also
reported; copy counters on Serial/OpenMP do not represent GPU traffic.
A/M times include finite checks, linear time includes A/M, setup includes
transfers, and transfer time can overlap those scopes. Do not add phase times.

## Qualification and reproducibility

`validation/kokkos_acceptance.json` declares controls and performance grids.
`test_kokkos`, `test_hispid_kokkos` and `test_by_kokkos` exercise the actual
compiled execution space. Original CPU Krylov controls and archived default
full-state/trace checks remain separate. Complete binary unknown/coefficient,
physical-field, gradient, independent FD-constraint and charge gates are
retained. The known HiSpID raw-P/coefficient comparison failure is not waived
by close physical fields, a new execution space or faster timing.

`validation/benchmark_kokkos.py` runs sequential fresh workers with three
interleaved repeats, both systems and both methods, at40×80×16,
80×160×16 and128×256×28. It records immutable source/build/image hashes,
hardware/affinity, actual runtime concurrency, true inner stopping histories,
phase/work counts, process memory and full states. Its32GiB explicit budget
is identical across variants. The frozen25ca064 reference receives only a
declared maximum-budget guard extension; its equations, allocation estimate
and stopping arithmetic are unchanged. The original8GiB capacity rejection
at the largest tier remains recorded.

No new spin.99 or Gamma10 physical investigation will precede the measured
performance report. Those are separate configurations selected by the user;
seed inputs and measured horizon properties will be reported separately.

## Prepared physical and consumer path (not numerically qualified yet)

`validation/run_extreme_kokkos.py` invokes the common physical validation
machinery with explicit CUDA selection, fixed GMRES forcing, separate output
directories and one retained grid per invocation. It rejects an incomplete
288-worker matrix before loading a native library. Its compilation receipt
must name `mcp__codex_app__compile_latex_document` and bind the successfully
compiled report and complete performance JSON by SHA256. The current compiled
98-worker snapshot cannot unlock this stage. A fresh final receipt is required.
`check_target_seeds.py --extreme` supplies new separate .99/Gamma10 controls
with refined, beam-aligned independent ADM quadrature. No new target has run.

Every extreme-study export is explicitly `diagnostic`: physical constraints,
refined charges, solved covariance and AthenaK horizons remain distinct gates.
The diagnostic-investigation switch records failed prerequisite gates and
does not change their criteria or inherit acceptance. The boost study retains
changed-separation configurations individually when calibrating d/measured
component Mirr. Its Fourier control uses a separate label.
Prerequisite bytes and all measured images are frozen and rechecked around
the solve/export. Exclusive case and row receipts bind the declared grid
prefix, controls and both raw NPZ artifacts; progress writes cannot erase
that binding. GPU UUID is retained per worker so an equivalent later
allocation can continue without rewriting provenance. Explicit execution
requires native convergence and weighted Linf within the declared tolerance;
an API success code alone cannot qualify a retained iterate.

Use `validation/portable_sampler_migration.py` with the checkpoint, exact
producer library, pure reference CPU consumer library, and a fresh output
filename. It compares identical coefficients in separate processes, checks
the exact continuous basis/maps, physical fields and metric derivatives at
off-grid and trial-horizon points, and binds artifacts and dependency images.
This is a sampler migration, without a PDE-equivalence or physical-acceptance
claim. A consumer with a mapped puncture Kokkos runtime is rejected, avoiding
a second Kokkos runtime in the AthenaK process. Actual numerical proof for each
new checkpoint is still required.

The charge and covariance replay tools accept separate `--results` and
`--raw-directory` inputs. New charge qualification uses an explicit bound
at most1e-5, at least five increasing radii, and separate polar-only and
azimuthal-only angular refinements. It checks native versus independent
integrals, fixed-origin transformation, angular changes in both finite-radius
values and intercepts, consecutive four-radius fit windows, and the reported
all-radius intercept against the finest window. Subsequent joint refinements
cannot bypass these checks. Qualified charge evidence remains separate from
binary acceptance.

For a new covariance study, provide a fresh `--output` and one
`--charge-evidence` file for each selected full grid. The physical free data
and sample-point mode must be identical and all three grid dimensions must
refine monotonically. Each fresh rotated solve repeats the bound native and
independent charge sequence at both displaced and global origins. It requires
those refinements as well as the existing covariance thresholds; the historical
12×24 preliminary integral cannot qualify this path. Explicit study stopping
checks require native convergence and the stated weighted tolerance. Rotated
and translated unknowns, diagnostics and histories are retained immediately
after solve, before measurement callbacks. Source files, arrays and images
remain bound by hashes throughout. The unchanged default invocation retains
its historical behavior. These new paths have metadata/refinement regression
tests but no fresh extreme numerical evidence yet.

AthenaK's binary driver accepts `--migration-proof` and an explicit
`--domain-half-width`. It rechecks checkpoint, executable, generated input,
proof and dependency hashes around every horizon worker. A separate `--common`
search also requires `--common-radius`; its configured center is recorded and
both modified balls must satisfy continuous retained-surface enclosure bounds.
Component searches keep their own masses/spins. All runs stay at time/cycle
zero. Failed common searches do not demonstrate absence.
