# HiSpID milestone status

Current status (2026-10-02): **exact single-hole targets pass; the moderate binary passes
preliminary constraints, charges and coordinate covariance**. Rest spin chi=.95 and lab speed v=.885 pass separately
and combined with generic directions. AthenaK's current isolated consumer
passes their initial-time FastFlow and three-level mesh constraint controls.
No evolution steps are taken. Solved-binary horizons and attenuation enclosure
have not yet been checked. Stronger accuracy and high-parameter binary
validation remain incomplete.

The previous native checkpoint is45bda9d. Default9cbf1108… uses regular
modal P and exact modal FD block preconditioning. The current experimental
charge/memory-only build126300dc… passes all native controls and17 Python
tests; separate-process fields, residuals and JVPs are bit-identical to9cb.
The default residual norm is unchanged. The sibling AthenaK consumer remains
bound to9cb; its current replay evidence is docs/hispid-current-controls.json.

The matched80×80→80×160×28 polar refinement improves fixed bulk physical
momentum RMS5.79x (.002435→.000421) and identical-grid g1 momentum RMS5.51x,
with decreasing auxiliary angular tails. Both binaries still fail acceptance.
The fresh80×160/104×208/128×256×28 sequence improves monotonically and meets
the preliminary finest-grid physical RMS limits. At128×256×28 nearH/M are
1.95e-6/6.02e-6 and bulkH/M1.25e-7/8.01e-5. The solve takes209.0s/39Krylov
iterations, process peak7.83GB. Refined angular charges now pass1e-7 at
allthree levels; finest extrapolated EPJ change is2.33e-7. The preliminary
physical/charge gate passes. Fresh three-grid solved coordinate covariance
also passes: finest metric/K rotation errors3.92e-11/2.84e-6, translation
errors<3.5e-13, fixed-origin charge error1.99e-4. Stronger accuracy and binary
horizon/enclosure checks remain pending.
See validation/polar_sequence_plan.json. Every failed result remains retained.

ADM integration has been sped up12.20x using tested analytic metric gradients.
The sphere grid now follows the prolate axis: both polar/azimuthal quadrature
and independent physical-FD flux comparisons pass their integration controls.
Earlier apparent polar convergence at fixed global-z phi64 was misleading:
phi128 shifted angular momentum by.00935. Those records remain unaccepted.
The integration repair alone does not promote a binary convergence gate.
The subsequent explicit three-grid integration check qualifies preliminary
moderate data only (`polar_sequence_refined_charges.json`). Ring reuse makes
64×64 extraction another21.12x faster, with differences<5.6e-13 and unchanged
ordinary sampler/operators. Its all-mode, rotated-frame, invalidation and
independent centered/off-center FD controls pass.

The optional sin6/(1-t)^6 row equilibration fails its declared exact-seed
far-source floor1e-14 (worst1.89e-13). No binary solve has used that norm.
The current default suite passes; the separately retained strong-anisotropy
restart32 inverse audit has7 failed rows. Restart128 passes its two difficult
m0 vector rows, including far-value/gradient checks, without changing the
original failed flags. Focused-map and wide-core experiments remain failed.

See docs/HISPID_PERFORMANCE.md for measured speedups and the integration
handoff; docs/VALIDATION.md and validation JSON retain full configurations,
source hashes, failure diagnostics and numerical history.

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
  the finest polar-refined sequence passes preliminary physical RMS limits;
  refined charge qualification and current solved coordinate covariance also pass.
  Earlier unregularized no-far-filter
  preliminary passes are historical and cannot satisfy the current gate.
- D: preliminary current off-grid convergence, refined centered charges and
  solved rotation/translation covariance pass. Stronger physical accuracy,
  refined fixed-origin charges and solved-binary horizon enclosure are pending.
  Historical checks retain their source fingerprints.
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
  combined run fails and is retained. Boost-only fixed-lmax
  quadrature to ntheta74 also passes: relative area
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
