# HiSpID serial performance and integration handoff

The solver and the AthenaK importer are implemented in isolated branches.
Exact isolated chi=.95 and v=.885 controls, separately and combined, pass
independent physical checks and AthenaK initial-time horizon controls.
The current regular-basis moderate binary passes preliminary physical and
refined charge, solved coordinate-covariance and direct AthenaK horizon/enclosure
gates. Stronger accuracy and revised high-parameter binaries remain pending. Its high-spin/boost range is
unvalidated; do not use these binary outputs as validated production data.

## Measured improvements

### Shared stopping norms and opt-in BY modal inverse (2026-10-02)

A new opt-in Fourier modal/block-tridiagonal BY preconditioner cuts the median
40×80×16 moderate-binary solve from **26.780s to 4.579s (5.85×)** while
retaining the original outer and inner stopping rules. Process peak RSS rises
from **91.66MB to 102.48MB (11.8%)**. Default BY and HiSpID states remain
bit-identical to the preceding archived production builds. Modal BY changes
its Krylov path but converges numerically to the same BY data.

The original criteria are not equivalent. BY stops on
`max |sin³(alpha) sin³(beta) F_raw|` and its linear absolute L2 target is
`1e-3 * Newton_max_residual`. HiSpID stops on the sixth-power version and
uses adaptive RHS-relative GMRES forcing. The old BY linear target in this
case is only about 4.8e-5–6.2e-5 of the RHS L2 norm, substantially tighter
than a fixed relative target of1e-3. Changing the outer weight alone also
changes HiSpID's adaptive forcing because it is computed from that norm.

The controlled comparison uses the SAME computational formula:

- Outer: `max_component max_node |(sin(alpha)sin(beta))³ F_raw,c| <=1e-12`.
- Inner: TRUE spectral `||F-J step||₂ / ||F||₂ <=1e-3` at every accepted step.
- Identical counts40×80×16, bare masses .6/.4, centers±3, generic lower spins
  and velocities, zero guess, serial CPU, restart64 for Hi GMRES and the
  inherited cap100 for BY BiCGStab. Exact input vectors are saved in the JSON.

Two interleaved fresh-process runs per protocol give these medians. MB is
1,000,000 bytes; RSS is captured before snapshots/physical verification.

| Protocol | Newton | Krylov | Spectral JVP | M applies | Solve | Ready for sampling | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|
| Hi native sixth-power/adaptive | 4 | 43 GMRES | 51 | 43 | 3.339s | 4.600s | 338.94MB |
| Hi cubic/adaptive (outer only matched) | 5 | 50 GMRES | 60 | 50 | 3.926s | 5.188s | 342.21MB |
| Hi cubic/fixed1e-3 | 4 | 42 GMRES | 50 | 42 | 3.262s | 4.528s | 319.14MB |
| BY inherited rules/line sweeps | 4 | 36 BiCGStab | 73 | 69 | 26.780s | 26.918s | 91.66MB |
| BY cubic/fixed1e-3/line sweeps | 4 | 27 BiCGStab | 62 | 54 | 22.582s | 22.718s | 91.55MB |
| BY inherited rules/modal | 4 | 13 BiCGStab | 32 | 24 | 4.579s | 4.718s | 102.48MB |
| BY cubic/fixed1e-3/modal | 4 | 10 BiCGStab | 25 | 17 | 3.844s | 3.984s | 102.51MB |

Both fully matched solves need four Newton steps. Krylov iterations are not
one-for-one comparable: BiCGStab generally uses two preconditioner/JVP pairs
per iteration, whereas GMRES uses one. BY modal has25 total spectral JVPs
including initial/true-residual checks, compared with Hi50. The BY core solve
is now only1.18× the Hi core solve; setup through first sample is0.88× Hi.
Matching tolerances alone reduces inherited BY iterations36→27 and time
26.78→22.58s. The similar modal inverse then reduces matched time another5.87×.
Per-step counts are BY lines7/7/6/7, BY modal2/2/3/3, Hi9/11/11/11.
All true relative linear residuals meet1e-3 (largest BY modal9.569e-4,
Hi8.242e-4); histories and native/unweighted max/RMS are saved explicitly.

The modal inverse computes a Fourier average of the POINTWISE current
potential `(7/8) A_BY² / psi⁸`, keeps full meridional blocks, uses the inherited
FD Fourier eigenvalue, and folds boundary ghosts exactly as legacy `Index`
(with no azimuthal parity). Its analytic orthogonal-prolate five-point stencil
avoids probing every original FD matrix column and analytically cancels mixed
terms. Factors are equilibrated and rebuilt at every Newton call, including
bare-mass adjustment. The original spectral F and JVP are unchanged.
A separate constructor using the stored original JFD remains the independent
control. Both inverses are checked against a full physical-grid even cyclic
projection of `SetMatrix_JFD`, including reflected faces/corners, rectangular
grids, every Fourier partner/Nyquist and refreshed unknowns/masses. Largest
inverse solution error is3.342e-14; normalized matrix residual1.772e-13.

Measured phase medians explain the speedup: inherited FD setup takes5.26s;
analytic modal setup takes0.200s. Inherited line preconditioning performs
13,800 full relaxation sweeps; modal inversion performs24 applications and
no sweeps. With the original rules, spectral JVP time drops7.67→3.38s through
fewer Krylov steps. Under matched rules, spectral JVP still takes2.65 of3.84s;
modal setup/application cost0.199/0.033s. The remaining BY cost is mainly the
original spectral Jacobian path, not factorization of the elliptic inverse.
The phase clocks are benchmark-only wrappers; they do not change arithmetic.

For original-rule modal BY, largest scaled nodal V difference is5.885e-15,
ADM mass differs3.109e-14 absolute, and puncture ADM masses are identical.
Fixed-forcing modal changes V by1.493e-14 scaled, ADM by7.17e-14 absolute,
and puncture masses by<=1.23e-14. Recomputed weighted spectral residuals
remain<=1e-12 in every run, with zero linear/modal failures. Nodal V/u values, scalar coefficients, and
six-point physical gamma/K/lapse/psi samples pass the predeclared1e-10 bounds.
All retained arrays are finite. Native Cartesian derivative arrays near the
punctures amplify roundoff and changed Krylov-path differences: matched modal
u_d23 differs by2.375e-4 scaled and u_d33 by9.306e-4 absolute. Largest differences occur at distance9.47e-7 from a puncture. Those
arrays are retained diagnostics, outside the1e-10 core-value/physical-sampling
gate; complete bitwise-state equivalence is claimed only for inherited defaults.
Fresh12×18×8 controls against original standalone `ec563aeb` pass for generic
spin/momentum, target ADM masses, Brill–Lindquist and the spin95 BY INPUT.
Defaults remain bit-identical to original; target-mass modal converges within
adm_tol1e-10 with numerically equivalent bare/end masses. These are solver
behavior controls, not new high-spin binary physical validation.

A common computational norm does NOT equate physical error. The maps,
seed geometry and one-versus-four PDE components differ. Independent
fourth-order Cartesian physical constraints use identical six points and
steps .004/.002/.001, with g=1 enforced for EVERY stencil sample. At step.001,
matched BY modal Hmax/Mnorm-max are4.055e-8/1.424e-9; matched Hi gives
7.043e-6/6.606e-3. The BY stencil reaches a differentiation/spectral floor;
Hi's momentum residual is stable across stencil steps, so this coarse grid
still has a material physical truncation error. This timing case is not an
accepted binary accuracy comparison or replacement for resolution sequences.
The independently computed H/M changes between old and modal BY are within
the declared FD-floor bounds.

Enable the improvement after setting default parameters:

```c
TwoPunctures_params_set_Int("TP_preconditioner", 1);
```

Legacy defaults remain0, with original200-sweep preconditioning and original
absolute linear gate. Requested modal requires scalar vacuum BY and even
nphi>=4. Optional `TP_linear_relative=1`, `TP_linear_rtol=1e-3` select the
matched forcing convention; they are not needed for the5.85× native-rule gain.
Invalid options, factor failures and failed true linear gates fail closed:
Newton does not accept the step, target masses are not updated, and a retained
failure reports diagnostics status1. The modal context and inherited BY global
parameter API are serial/non-reentrant. No conditioning guarantee is inferred
from exact-zero pivot checks.

Hi default row power remains6. `HISPID_EXPERIMENT_FLAGS=-DHISPID_ROW_POWER=3`
builds the explicit cubic experiment (distinct residual-scaling tag).
`HiSpID_solve_with_forcing` / Python `solve(linear_rtol=.001)` supply fixed
forcing without changing Config/Diagnostics layouts. New work/history queries
record adaptive or fixed solves; buffer capacity is validated. The positive
row-scaling control finds residual/JVP differences1.013e-17/8.645e-19 scaled,
with physical fields bit-identical between powers3/6 for IDENTICAL INJECTED
unknowns (not separate solved states), and all default operators
bit-identical to the archived sixth-power producer.

All native controls,127 modal/invalid/failure checks,2056 line-cache checks,
24 GMRES/history controls,30 Python tests and physical API/lifecycle regression
pass. The sole independent formulation reviewer checked the derivation,
boundaries, signs, guards and validation gates. Evidence:
`validation/by_modal_acceptance.json`, `validation/common_stopping_moderate_40.json`,
`validation/stopping_row_scaling_controls.json`, and
`validation/by_modal_verification_v2.json` and
`validation/modal_final_fingerprints.json`. Raw states, logs, frozen precision
sources, object/image fingerprints and failed launch diagnostics are retained
locally. Prior physical producer hashes and all failed high-parameter gates
remain unchanged; no AthenaK consumer/evolution change is made here.

Reproduce in fresh absolute build directories, with one CPU worker at a time:

```sh
TPROOT="$PWD"
MODAL_BUILD="$TPROOT/build-by-modal-replay"
make -j1 OBJD="$MODAL_BUILD/obj" LIBD="$MODAL_BUILD/lib" \
  HISPID_DIR="$MODAL_BUILD/hispid6" test-hispid test \
  TEST_EXE="$MODAL_BUILD/test_physical_api.x" "$MODAL_BUILD/lib/libTwoPunctures.so"
make -j1 OBJD="$MODAL_BUILD/obj" LIBD="$MODAL_BUILD/lib" \
  HISPID_DIR="$MODAL_BUILD/hispid3" HISPID_EXPERIMENT_FLAGS=-DHISPID_ROW_POWER=3 hispid
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONPATH="$TPROOT/python:$TPROOT/validation:$TPROOT/examples"
python3 validation/build_by_solver_timer.py --objects "$MODAL_BUILD/obj" \
  --output-dir "$MODAL_BUILD/timing"
python3 validation/benchmark_common_stopping.py \
  --hi6 "$MODAL_BUILD/hispid6/libHiSpID.so" --hi3 "$MODAL_BUILD/hispid3/libHiSpID.so" \
  --by-library "$MODAL_BUILD/timing/libTwoPuncturesTimed.so" \
  --output validation/common_stopping_replay.json
```

### Behavior-preserving memory and line-solve changes (2026-10-02)

Two alternating fresh-process measurements per version use the same lower-spin
unequal-mass configuration, 40×80×16 counts, zero guess, one CPU thread and
native tolerance 1e-12 on the Apple M5 Pro. The following are medians;
MB means 1,000,000 bytes and memory is whole-process peak RSS captured before
verification snapshots. Timing includes identical full-precision BY logging.

| Backend | Solve before → after | Peak RAM before → after | Result |
|---|---:|---:|---|
| HiSpID | 3.332s → 3.315s | 496.45MB → 338.35MB | 31.8% less RAM; speed essentially unchanged |
| Bowen–York | 57.338s → 26.988s | 86.28MB → 91.62MB | 2.12× faster; cache costs 6.2% more RAM |

Setup through first sample changes from 4.590s to 4.579s for HiSpID and
57.477s to 27.129s for BY (2.12× faster). Two repeats establish the observed
workload cost, not a scaling law or a statistically significant HiSpID speedup.

HiSpID allocates GMRES V/Z columns only when visited, keeps them across restarts,
and removes an unused full-size scratch vector. Restart64, iteration limits,
zero initialization, orthogonalization, update order and stopping rules are
unchanged. Savings depend on the largest visited Krylov bank; a solve using
all restart columns will retain nearly the previous history storage. The
conservative worst-case context memory guard is retained.

BY caches nonsymmetric tridiagonal factors and off-line entry masks for each
fixed finite-difference Jacobian. Each line owns its mutable RHS/solution
workspace. Factors are rebuilt for every new Newton/BiCGStab call, including
the target-mass loop. The implementation keeps GSL's division and elimination
expressions, stencil subtraction order, reflected endpoint assignments,
plane/color/line update order, all 200 relaxation sweeps, and the spectral
Newton/BiCGStab operations. Unsupported pivots, oversized row masks or failed
cache allocation select the inherited GSL path.

The baseline BY equations, coordinates, spectral operations, Newton and
BiCGStab source match original standalone TwoPunctures commit
`ec563aeb672235b9443c330f9cde65f7246e8ea4`. A fresh build of that commit gives
bit-identical v/u/cf_v arrays, stored and recomputed F, and Cartesian spectral
correction samples at the full benchmark grid. Four additional 12×18×8
controls cover generic spin/momentum, target ADM masses, a zero-source binary,
and the spin95 BY input. All reach the native tolerance and match the original
state exactly. Target-mass final bare masses and internal-end ADM masses also
match exactly. This checks the first fresh solve; the fork's existing explicit
ownership/lifecycle guards intentionally differ from upstream static reuse.

Both baseline and optimized full-grid BY runs have identical full-precision
Newton/BiCGStab traces, 37 snapshot arrays, ADM diagnostics and physical samples.
HiSpID's unknowns, recomputed residual, physical metric/curvature/gradient
samples and Newton/Krylov counts are bit-identical. Separate controls include
2,056 BY checks, comparing legacy-GSL/cached results after every step of 800 relaxation
sweeps on refreshed rectangular-grid Jacobians, and 18 eager/lazy GMRES cases.
The final complete native suite, 30 Python tests and physical API/sequential
lifecycle regression pass. The sole formulation reviewer found no defect.

BY retains its original max-absolute Newton residual weighted by
`sin^3(alpha) sin^3(beta)`, and the unnormalized L2 BiCGStab tolerance
`1e-3 * Newton_max_residual` with cap100. HiSpID uses
`sin^6(alpha) sin^6(beta)` with its GMRES forcing and damped Newton rules.
The numeric tolerance does not represent the same stopping norm or physical
accuracy between backends. Neither set of criteria was changed here.

The implementation is commit `a4fbd5f`. Evidence is
`validation/solver_efficiency_moderate_40.json` and the stronger
retained-record checks in `validation/solver_efficiency_verification.json`.
They record explicit library hashes, configuration/diagnostic equality,
array shape/dtype/bit checks, trace hashes, source hashes and original-code
controls. Raw states/logs are retained in the ignored validation/raw directory.
Old producer libraries/checkpoints and failed physical gates remain bound to
their original hashes; these performance tests confer no new binary physical
acceptance.

For a fresh serial build and all controls, use fresh output directories:

```sh
make -j1 test-hispid test \
  HISPID_DIR="$PWD/build-efficiency-replay" \
  OBJD="$PWD/build-efficiency-replay/obj" \
  LIBD="$PWD/build-efficiency-replay/lib" \
  TEST_EXE="$PWD/build-efficiency-replay/test_physical_api.x"
```

`validation/benchmark_solver_efficiency.py` runs one numerical worker at a time,
alternates before/after order and requires explicit frozen libraries. Supply
both a pre-change snapshot (commit c4158cb) and current libraries, a fresh
`--output`, and `--grid 40:80:16 --case moderate --tolerance 1e-12 --repeats 2`.
The timing wrappers replace TP_Newton.o and alter only printf precision, as
specified in the verification provenance. `validation/compare_by_upstream.py`
uses only original-compatible symbols. `validation/verify_solver_efficiency.py`
rechecks retained arrays, configurations, diagnostics and traces without
re-solving or overwriting evidence.

### Cold comparison with the local Bowen–York backend

The lower-spin unequal-mass case uses masses.6/.4, centers(±3,0,0),
dimensionless rest seed chi=.390512/.374166 and the generic spin/velocity vectors in
`examples/configs.py:moderate`, with far filtering disabled. Both backends
start from zero at40×80×16 and use one CPU thread on the Apple M5 Pro,
with-O3 builds and nominal native weighted-residual stopping tolerances.
BY inputs match the isolated seed lab charges:
`P=m Gamma v` and
`S_lab=Gamma S_rest-Gamma² v(v·S_rest)/(1+Gamma)`.
The nodes/maps and geometries differ; this is a comparison of grid counts
and bare/seed-charge inputs, not equal measured physical accuracy or horizon
properties. No binary acceptance is inferred.

| Nominal tolerance | HiSpID solve | BY solve | HiSpID total to first sample | BY total to first sample | Both reached tolerance? |
|---|---:|---:|---:|---:|---|
| 1e-10 | 2.606s | 51.056s | 3.913s | 51.200s | Yes |
| 1e-12 | 3.359s | 56.833s | 4.633s | 56.972s | Yes |
| 1e-14 | 3.661s | 202.666s | 4.964s | 202.819s | No: BY residual2.37e-13 |

The archived pre-optimization1e-12 pair makes HiSpID16.9× faster for the elliptic solve
and12.3× faster through setup and first sampling. Peak process RSS is
495.94MB versus85.54MB, about5.8× higher. The1e-10 row uses qualified
producer126300dc; the1e-12/1e-14 rows use the later exact-reuse sourceb96ee4b1.
No accepted performance ratio is assigned to the failed1e-14 pair. Each row
is one paired measurement, not repeated-run uncertainty or a scaling curve.
The80×160×28 BY screen was interrupted; its completed HiSpID time41.51s
is retained separately and cannot provide a completed ratio.

BY solves one scalar equation with analytic momentum, while HiSpID solves
four coupled curved equations. BY uses its unchanged200-relaxation-sweep
preconditioner; HiSpID uses exact modal FD block elimination. The ratio
applies to these local serial implementations. It does not predict the cost
of another BY implementation or finer grids.

`validation/by_comparison_summary.json` binds the full configurations,
loaded libraries, timings, memory, raw-file hashes and limitations. The final
1e-12 pair uses spectral first sampling for both. The first1e-10 and failed
1e-14 rows retain BY's default Taylor first sampling; this distinction does
not affect solve times. BY eagerly constructs coefficients inside its
initial-data call; HiSpID's lazy coefficient construction is included by
the first-sample metric. Verification, ADM extraction and horizon searches
are excluded. All failed and interrupted records remain separate.

For a fresh measurement of the current source, link the timing wrapper instead of the legacy Newton object;
it includes the original Newton source without changing its mathematics:

```sh
mkdir -p build-hispid-by-replay
make -j1 "$PWD/lib/libTwoPunctures.a"
cc -std=c99 -O3 -fPIC -Iinclude $(gsl-config --cflags) \
  -c validation/by_newton_timer.c -o build-hispid-by-replay/timer.o
cc -shared -O3 -fPIC build-hispid-by-replay/timer.o \
  obj/TwoPunctures.o obj/TP_CoordTransf.o obj/TP_Equations.o \
  obj/TP_FuncAndJacobian.o obj/TP_Utilities.o $(gsl-config --libs) \
  -o build-hispid-by-replay/libTwoPuncturesTimed.so
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  PYTHONPATH=python:validation:examples python3 validation/benchmark_bowen_york.py \
  --hispid-library "$PWD/build-hispid-preconditioner-reuse/libHiSpID.so" \
  --by-library "$PWD/build-hispid-by-replay/libTwoPuncturesTimed.so" \
  --case moderate --grid 40:80:16 --tolerance 1e-12 \
  --output validation/by_comparison_fresh_replay.json
```

Supply a fresh output label and the explicit isolated HiSpID build path.
Historical source fingerprints remain bound to their own measurements.

- At identical moderate80²×28 free data, exact block elimination of the
  five-point modal FD preconditioner reduced the solve from296.95s/896
  Krylov iterations to31.33s/64 iterations:9.48x. The full residual and JVP
  remain pseudospectral; this is a preconditioner improvement. Physical
  failure flags are retained (`modal_block_equivalence.json`, `results.json`).
- ADM energy now uses the independently tested analytic physical metric
  gradient. At matched16×32 sphere nodes and R40/200/1000, integration time
  fell13.657s→1.120s (12.20x). Energy differences are<2.7e-11; P/J are
  bit-identical when grid points match (`analytic_charge_method_comparison.json`).
- Aligning the sphere polar axis with the prolate axis removes the severe
  azimuthal alias of the retained meridional polynomial. On the same failed
  80×160×28 binary at R200,256×64→256×128 changes EPJ by6.19e-14;
  256×128→384×128 changes it by6.35e-10. Independent Cartesian-FD fluxes
  and centered tensor covariance pass1e-8 at R40/200/1000; maximum
  discrepancies are5.52e-11 and1.61e-10. This validates integration, not
  the binary geometry (`aligned_charge_fd_comparison.json`,
  `charge_quadrature_polar2_aligned.json`).
- The allocation guard now counts the current five-point stencil instead
  of the obsolete76-point stencil. The default2048MiB and explicit8192MiB
  maximum are unchanged. The new estimate also retains dense radial block
  factors, their shared Fourier/vector groups, and the GMRES basis.
- Centered aligned sphere rings now reuse P/Pt/Peta across phi. Matched
  64×64 quadrature at R40/200/1000 on128×256×28 data falls22.081s→1.046s:
  another21.12x, with maximum EPJ difference5.52e-13. Ordinary sample,
  residual and JVP paths remain bit-identical. Phase-shifted all-mode
  first-gradient, translated/rotated-frame and changed-unknown controls
  pass; independent centered/off-center FD fluxes agree within5.82e-11.
  Off-center spheres retain ordinary interpolation.

A subsequent exact serial reuse removes repeated azimuthal averaging per
Fourier row and moves the vector LU/transfer bank between Newton steps.
Scalar potential/factors and every row scale still rebuild. On a matched
fresh rotated80×160×28 solve, time falls59.97s→55.44s (1.082x in one paired
measurement), with identical4Newton/45Krylov counts and bit-identical final
unknowns, physical fields and residuals. Peak RSS is recorded per worker.
The rectangular private control proves bitwise matrix/factor/random-RHS
equivalence, changing scalar potential, first/second factor counts7+7/7+0,
malformed-cache rejection and direct modal inverse error5.9e-18. A scaled
Fourier round-trip fixture failed identically in cached/uncached paths at
6.5e-10; that first diagnostic remains retained. Separate-process field
checks on allthree saved grids and80×160 residual/JVP checks are bitwise.
All native and17Python controls pass. These source changes do not relabel
the original binary acceptance records.

On the current80×160/104×208/128×256×28 moderate sequence, all near/bulk
physical H/M RMS values improve monotonically. The finest nearH/M are
1.95e-6/6.02e-6 and bulkH/M1.25e-7/8.01e-5; the209.0s solve uses39Krylov
iterations on one thread. Each grid's ADM extraction uses polar orders
2Npolar/3Npolar and phi64/128 at R100/200/400. All satisfy the declared
1e-7 integration bound. The two finest extrapolated EPJ vectors differ by
2.33e-7. `polar_sequence_refined_charges.json` records this preliminary
qualification and its explicit failed stronger/full-binary flags.

Fresh rotated/translated solves of the same sequence pass the declared
coordinate gate. Finest metric/K rotation errors are3.92e-11/2.84e-6;
translation errors are<3.5e-13. Centered charge error is3.40e-10 and fixed
global-origin charge error1.99e-4 versus the preliminary.005 bound. The
fixed-origin extraction still uses coarse12×24 integration and needs
separate angular refinement for stronger claims. The whole covariance
process peaks at10.32GB over sequential contexts (one CPU thread), while
the original finest solve alone peaks at7.83GB. Context budgets do not
bound allocator retention or whole-process RSS.

Jobs and builds use one CPU thread. No OpenMP/MPI/GPU implementation or
scaling claim is made. CPU time and peak memory must justify that next
step. Perlmutter SSH access is available; no allocation has been made.

## Review and replay

Native branch `codex/hispid` lives in the sibling `TwoPuncturesC` worktree,
based on68287742f4920f4ea39b7dac1571c81eefe2ff8f. BY equations remain
unchanged; TP_Newton now caches fixed-JFD line factors as tested above.
`docs/HISPID.md` contains build, ABI, conventions and replay instructions;
`docs/VALIDATION.md` and retained JSON give configurations and failed gates.
Always supply an absolute native-library path. A checkpoint is bound to its
source SHA, basis identifier, maps and full grid. Normal replay rejects a
mismatch. Explicit separate-process bitwise witnesses may authorize an
API-only initial guess or read-only integration check; fresh solves and
acceptance checks are still required. They never relabel a failed binary.

AthenaK branch `codex/hispid-pgen` lives in the sibling `AthenaK` worktree,
based exactly on PR790 head22baa243970fa1880b2bbc48e88a590069d55e47 from
`HengruiZhu99/athenak:project/z4c_overhaul`. Its `z4c/hispid` pgen loads
portable physical gamma/K data, fills active and ghost cells, and checks
ADM/Z4c round trips. `docs/hispid.md` gives its build/replay commands;
`docs/hispid-current-controls.json` binds the current executable and exact
seed tests to native9cbf1108… source provenance. Direct native-geometry
finder controls and independently refined mesh-constraint checks are
reported separately. These controls take zero evolution steps.

The checked128×256×28 moderate checkpoint has now been exported and tested
with the separate exact-source AthenaK consumer. Both component surfaces pass
expansion RMS1e-7 at lmax16,ntheta32/48. Fixed-order area changes are<7e-12
relative. Real-harmonic coefficient bounds certify complete inner g/operator
balls on the retained surfaces, with margins.176225M/.115817M after an
empirical refinement allowance. Upper-radius bounds certify distinct
components. The mesh import round-trip error is4.1e-16; no evolution steps
are taken. See validation/polar_sequence_horizons.json and the sibling
AthenaK docs/hispid-moderate-binary.json for source/consumer fingerprints.

This direct-geometry finder check is separate from mesh-resolved finder
accuracy. The noncompact f/F attenuation tails require exterior constraint
checks and cannot be described as entirely inside a finite horizon. Reported
coordinate spin is not a generic approximate-Killing-vector spin. Stronger
physical accuracy and revised high-spin/boost binaries remain unvalidated.

Both isolated branches were pushed to the user's forks on2026-10-02.
No merges, shared installations, main-project branch changes or production
evolutions have been performed. GPU work remains last priority.

The revised rest-spin chi=.95 binary now has three successful serial Newton
solves. Its128×256×24 solve takes219.32s with4Newton/54Krylov iterations;
the whole solve/verifier process peaks at6.74GB. Physical constraints still
fail, so timing is a diagnostic workload rather than an accepted high-spin
benchmark. The exact126300dc producer exports the labeled checkpoint and
AthenaK imports it with ADM/Z4c round-trip error4.11e-16.

Both coarse horizons are measured, with each Christodoulou mass.5006533,
area8.3135384 and coordinate chi.9475223. A measured-radius initial guess
reduces the matched coarse search213.69s→44.39s (4.81x), with relative area
change9.4e-13; this changes the initial guess only. The complete bounded
measurement/failed-fine-control process takes454.97s and peaks234.11MB on
one CPU. Fine lmax16 expansion stalls above the strict threshold, with
angular variance near2.32e-7. Fixed-order quadrature/higher-order checks
remain pending. See validation/spin95_local_horizons.json. No new parallel
or accelerator work is justified by these timings alone.

The equation/Krylov separation and complete 2×2 binary comparison are recorded
in [KRYLOV_MATRIX_RESULTS.md](KRYLOV_MATRIX_RESULTS.md). Both defaults are
bitwise preserved. Hi BiCGStab reduces RSS but fails strict saved auxiliary-field
and coefficient equivalence; physical agreement is reported separately and the
backend remains opt-in.
