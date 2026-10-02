# HiSpID milestone status

Current development status (2026-10-02): the committed sampler/API build
`2dbf450d…` has accepted exact seed and AthenaK initial horizon controls at
chi=.95 and v=.885 (including their combination). Binary acceptance with
uniform axis regularity remains pending. The current uncommitted native build
`9cbf1108…` uses regular modal P on analytically mapped Gauss grids, with a
direct modal FD block preconditioner. Seventeen Python controls and independent
all-jet AD controls pass. Fresh moderate24²×12/40²×20/56²×28 solves converge
internally but **fail** independent physical/charge gates. At56²×28 near
H RMS2.79e-5,M RMS5.32e-4; bulkH6.35e-5,M1.63e-3. The first80²×28 run
failed at the Krylov limit after177.7s. The cache/restart80 retry converged
internally in296.9s but failed the physical gate. Exact modal block elimination
reduced this to31.3s and64 Krylov iterations (previously896), retaining the
same failed physical residuals. At104²×28 the solve takes42.4s/54 Krylov
iterations; nearH/M RMS5.72e-6/9.49e-5 and bulkH/M3.98e-6/2.31e-4 still fail.
Unfiltered152²×28 reaches nearH/M9.76e-7/1.66e-5 and bulkH/M4.80e-7/1.74e-4;
far40 has essentially the same failed momentum norms. The current build adds
constant-annihilating meridional row differences. All17 Python controls and
the strict analytic flat-source inversion audit through64² pass. The prior
64² constant-source inverse missed its strict linear tolerance while its
reconstructed field/gradient error was8e-18; that outcome is preserved.
Off-axis sampled fields remain bitwise identical; residual/JVP changes are
5.81e-16/2.43e-17 in separate-process replay. Fresh current-basis56²/80²
solves also fail:80² nearH/M RMS2.13e-5/2.48e-4 and bulkH/M
2.87e-5/2.44e-3, matching the earlier direct-block result. These rounding
improvements do not promote prior failed binaries.
The cache-reordered derivative application and centered
preconditioner projection preserve off-axis fields, residuals and JVPs
bit-for-bit in isolated old/new library comparisons. Native controls now
include random-mode sampler/Cartesian-gradient consistency and zero
constant-mode leakage through the preconditioner. Random-vector modal FD
inversion through64² agrees to2.5e-14. Independent Cartesian scalar/vector
manufactured checks pass all axis segments and ordinary transverse points.
No current regular-basis
binary is accepted.
The default full native/Python suite passes. A read-only map-query ABI and
exact experimental-map identifiers leave default-map fields, residuals and
JVPs bitwise identical in separate processes (`collocation_map_api_equivalence.json`).
An isolated build-map experiment, radial stretch.05 and angular stretch3,
is retained as failed. Its source is `e7a88824…`; all three80²/104²/128²
solves miss the requested1e-14 internal tolerance at line search and fail
physical acceptance. BulkM RMS is.00110/.0101/.00261. Its basis token
encodes the exact constants; checkpoints
cannot be replayed across maps. It targets resolution of the narrow g/f
cores and is not an accepted result. Map changes require fresh derivative,
physical, charge-quadrature and convergence checks. All17 Python controls and
finest64² Cartesian scalar/all-three-vector operator controls pass. Strict
scalar N64 constant-mode relative inversion audits fail at3.90e-13 and
5.98e-13 versus1e-13, with physical correction value/gradient errors<2.2e-16;
failed flags remain unchanged. Random degree39 all-Cartesian-Hessian checks
pass at1.75e-13 normalized discrepancy. Coupled flat Navier inversions
through32² pass their declared1e-11 relative/field bounds; maximum sampled
field/gradient discrepancy is1.65e-14.

Replaying the same N80 retained polynomial on104²×28 without a solve exposes
g1 momentum component RMS up to4.46e-4, independently of the Cartesian FD
verifier. Doubling only phi to56 leaves RMS<5e-12. This localizes the gap
to continuous meridional PDE representation; true original-node unweighted
residuals are small and do not support large ignored core equations as its
sole cause. Separate polar-only/radial-only doubling gives g1 momentum
component RMS up to1.29e-4/3.07e-5, so both directions contribute.
The separately labeled wideg/actual-operator N56 run fails at the Krylov
limit; N80 was deliberately interrupted and N104 was not attempted. The
matched flattened-core N80/N104 runs also fail physical acceptance, with
bulk momentum RMS3.58e-4/4.16e-4. Their g-ball screen uses.2/.8 of each
contracted isolated Kerr throat and still needs computed binary horizon
enclosure. No configuration is accepted by these changes.
A constant-anisotropic-metric inverse audit (no black-hole sources) passes
all36 weak-anisotropy rows, but7 of36 stronger-anisotropy rows fail with
restart32/max300 Krylov. Its exact-P action agrees with an independent
Cartesian RHS to roundoff. A separate restart128/max1024 diagnostic passes
the two difficult m0 vector rows in82/85 iterations, with sampled errors
1.47e-12/5.56e-15. This confirms finite-restart stagnation in that control;
no linear tolerance or original failed flag is relaxed.
Previous binary results below are explicitly historical (unregularized basis).

Baseline: `68287742f4920f4ea39b7dac1571c81eefe2ff8f`. Branch: `codex/hispid`.
Worktree: `/Users/hz0693/research/lazarus/.hispid-worktrees/TwoPuncturesC`.
Original native checkout was clean; Python `main` initially had staged work.
Both remain untouched by this task while the main project continues its work.
No shared AthenaK files, installations, or simulations will be changed. Jobs use
one CPU thread; builds use `make -j1`. No merges into the main project.

- Isolation and baseline inspection: complete.
- Reference formulation and original implementation investigation: complete;
  local thesis access resolved. Historical stuffing and spin conventions differ.
- A: geometry/operator implementation and manufactured validation: passed.
- B: independent single-hole and boosted seed validation: passed.
- C: current regular-basis moderate unequal-mass, generic-spin/boost binary:
  independent acceptance remains failed. Earlier unregularized no-far-filter
  preliminary passes are historical and cannot satisfy the current gate.
- D: current independent off-grid convergence/axis/charge/covariance gates:
  pending. Historical solved rotation/translation checks passed; they do not
  accept the current basis. Refined charge angular quadrature is required.
- E: HS99UU160²×24 passes the finest-grid local strict physical thresholds
  and the declared energy comparison. Its original all-norms monotonic gate
  remains failed because bulk H is below the calibrated verifier resolution.
  Direct Gamma=sqrt5 binaries and flattened-operator continuation fail;
  a separately labeled actual-metric operator control reaches the target
  Newton tolerance but fails physical acceptance: at112²×8 near M RMS is
  1.35e-3 and the final ADM-energy grid difference is.025M.
  Exact published high-boost reproduction remains incomplete: bare masses,
  companion widths and original parameter files were not recovered.
- Revised user targets: rest seed chi=.95 and lab speed=.885, including a
  generic combined binary. Independent binary acceptance is still pending.
- AthenaK: isolated branch `codex/hispid-pgen` starts at PR790 head
  `22baa243970fa1880b2bbc48e88a590069d55e47`. A new physical-data pgen builds,
  imports portable checkpoints and checks ADM/Z4c round trips. It can call
  FastFlow at cycle/time zero and optionally compute mesh constraints.
- Exact horizon controls: Schwarzschild from two guesses and stationary
  Kerr chi=.95 pass; Kerr area/spin agree to roundoff. Stable harmonics let
  boosted Schwarzschild v=.885 pass at lmax48,ntheta50: expansion RMS9.99e-8,
  sampled shape relative error3.78e-8. Coarse angular strict flags remain false.
  Combined chi=.95,v=.885 seed passes with flow alpha=.2 at lmax48,ntheta50:
  expansion RMS9.92e-8, area relative error1.75e-14 and sampled shape relative
  error3.35e-7. Fixed-lmax quadrature refinement to ntheta74 changes its area
  by8.44e-15 relative and expansion RMS by3.61e-12. The original alpha1
  combined run fails and is retained. Boost-only quadrature and binary
  boost-only fixed-lmax quadrature to ntheta74 also passes: relative area
  change1.44e-15, RMS change1.79e-12. Binary enclosure remains pending.
  These controls use exact isolated seed data. Independent generic target
  seeds pass near/bulk constraint RMS<1.4e-8 and charge error<2.3e-7 after
  extraction-radius refinement to10240; short-radius failures are retained.
- Horizon verification of binary attenuation: not yet established; no
  exterior-vacuum claim for a binary is warranted.

## Acceptance gates (declared before binary/high-parameter runs)

Use geometric units and normalize lengths by total seed rest mass M. Test
physical H=R+K^2-K_ij K^ij and physical contravariant M^i=D_j(K^ij-gamma^ij K)
with an independent Cartesian finite-difference implementation, never the solver
equations alone. Record rms and maximum separately in the g<1 region, within
2M of a hole outside that region, and farther exterior.

1. Analytic/manufactured derivatives, curvature, and longitudinal operators:
   errors below 1e-8 at ordinary points, with finite-difference refinement.
2. Exact seed controls: exterior physical rms H and M components below 1e-7
   in units M^-2; independent derivative refinement must reduce errors until
   cancellation dominates. Schwarzschild/Kerr ADM E,P,J extrapolated over
   increasing radii must agree to 2e-4 relative or 2e-5 absolute.
3. Moderate binaries: scalar and momentum internal residual each below 1e-9
   after row scaling; positive spatial metric and conformal factor; independent
   off-grid near/exterior rms constraints below 1e-4, with at least three
   spectral resolutions showing improvement and stable charges (<0.5%). This
   is an initial validation gate, not the published precision of 1e-9--1e-12.
4. Rotation/translation covariance: seed tensors to 1e-11 relative, solved
   fields/charges to measured spectral truncation error; translated angular
   momentum follows J'=J+c cross P.
5. Only after 1--4 pass, run higher-spin/boost cases and characterize failure
   boundaries. Published reproduction requires specified gauge/attenuation,
   convergent independent constraints, and matching charge tolerance; horizon
   mass/spin comparisons remain unverified without an apparent-horizon solver.

Before the first high-regime run, the stronger physical threshold is fixed at
near and bulk RMS 1e-6, with maxima below 1e-4 and improvement on three fine
grids. Adaptive Cartesian verifier steps are 0.001 times distance to the nearest
puncture, capped at 0.002M, and checked at factors 0.5 and 2. HS99UU uses the
thesis's stated bare inputs, QI conformal factor, f width 0.2, and no g or F.
Its ADM energy target is 0.980124, with a comparison tolerance of 2e-4M; this
comparison cannot establish exact historical reproduction because the modern
trace projection differs. The high-boost local benchmark targets convergence
and charge stability only: no exact published charge target is available for
its full bare configuration. Angular quadrature and radial extrapolation must
be separately checked before accepting a charge comparison at this precision.

Failed-case configurations, iterates, diagnostics, and runtimes are retained.

## Historical evidence (unregularized basis)

The initial far40 grids 12²×8,20²×12,28²×16 failed physical validation and
are preserved in `validation/failed_moderate_coarse.json`. Independent review
found and verified a repair for the off-grid cosine Nyquist derivative.
An exact Brill–Lindquist far-filter control exposed underresolution of the
known scalar shell; the exact variable change `u=W+v`,
`W=sum((1-F)(psi_seed-1))`, now carries that term analytically without changing
the free data or constraints. Its exact-control and analytic-JVP tests pass.

No-far-filter near H RMS decreases from 3.25e-4 to 1.20e-4 to 1.50e-5;
near momentum norm RMS decreases from 2.62e-4 to 7.13e-5 to 9.12e-6.
At 56²×28 bulk RMS H=1.14e-8 and M=6.77e-8; extrapolated ADM energy is
.9944889342. This passes the declared 1e-4 preliminary gate, not the stronger
1e-6 near-hole target or publication-level accuracy. No high-regime case has
passed the original aggregate gate. The stable seed-divergence identity
repaired the first HS99UU source's momentum cancellation. Angular refinement
160²×16 to160²×24 reduces near H RMS from1.55e-6 to3.18e-7 and M RMS from
8.69e-8 to2.70e-9. Exact vacuum controls at the same18bulk points match
the observed bulk H stencil noise; raw gates remain unchanged.
See `validation/verifier_floor_highspin.json`.
Refinement explicitly uses4096 or6144 MiB context budgets on the48 GiB host;
the default remains 2048 MiB and only one numerical job runs at once.
Large global raw conformal extrema close to singular map
coordinates are retained; physical-equivalent collocation constraints at
g=1 are now reported separately. They do not replace independent checks.

The sampler/API build is SHA256
`2dbf450d83d87ad58484fe9b676327dedd93ae8d185759281898322c338b49cb`.
The archived pre-axis build is79e96b4c…; off-axis physical fields and the
moderate residual/JVP compare bit-for-bit in fresh, separate processes in
`validation/axis_api_migration.json`. The original same-process comparison
was invalid: dyld coalesced the two images by their identical install name.
It is preserved with an invalidation note; affected revalidation results
were withdrawn and rerun. Python checks the actual symbol image with dladdr
and freezes a process-wide first-load SHA. Distinct builds must be compared
in separate processes. Fresh current-build moderate/covariance gates now
pass under `moderate_sampler_api`; old checkpoints are only an explicitly
proved initial guess. Sampler and
analytic physical metric-gradient controls pass, including fixed W and
rotated frames. A finite spectral expansion can violate map-axis endpoint
regularity: solved-axis Hamiltonian tests can scale as h^-2. The manufactured
axis control passes, but does not remove this solved-data limitation. This
needs resolution characterization or a regular basis before a uniform
continuum accuracy claim. Point-only sampling changes do not repair it.
