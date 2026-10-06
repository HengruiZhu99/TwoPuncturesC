# Trumpet implementation and acceptance ledger

Goal started 2026-10-05. Native base4ce5974; AthenaK consumer base9495fdb9.
Dedicated branches: codex/hispid-trumpet and codex/hispid-trumpet-pgen.
The prior worktrees and raw results are unchanged. No previous benchmark
campaign is resumed. Current status: seed-family API, context/cache/sampler
dispatch and versioned checkpoint integration implemented; CUDA controls and
full consumer build pass. Isolated spin.99 and the finest moderate binary
horizon checks pass, including angular convergence and modified-region
enclosure. Two grids in the declared sequence meet the physical error bounds,
but the sequence fails decreasing bulk errors. No full binary acceptance yet.

## Pinned geometry and conventions

AthenaK project/kerr-trumpet-spin09 at81c7d8030c49489577c54b977b8b7f7b42d9a96e;
provider blob f297008baa186b98ee6899846365e2dba7102e73, also present in
project/tde atd2ef71e9797f194a9631d78f8fde885638c28b19.
Source: Dennison, Baumgarte and Montero, arXiv:1409.1887v2, Eq.12 withR0=M.
The reference has positive stationary ADM lapse, not ordinary stationary
1+log slicing. Residual evolution background subtraction is not part of
the elliptic constraint equations.

Use signature(-+++), Kij=-Lie_n(gammaij)/2 and rest spinS with a=|S|/m.
LetR=r+m, Sigma=R^2+a^2 cos(theta)^2 and
X=(R^2+a^2)^2-a^2 r^2 sin(theta)^2. The existing Cartesian expressions
generalize by replacing the z axis with the unit spin vector. For any boost,
pull back the complete four-metric by the Lorentz map and differentiate its
time-dependent spatial block. Derive K from the complete ADM identity.
Retain the q^2-alpha0^2 v_i gamma0^ij v_j>0 spacelike check and positive
spatial principal minors; never clip a failed slice. Investigate an alternative
slicing only if this construction fails the required domain.

The CTT scalar choices are(Sigma/r^2)^(1/4) or det(gamma)^(1/12).
Both have leadingr^-1/2 behavior with directional coefficients. CTT covariant
A=psi^2(Kij-gammaij K/3), while the AthenaK Z4c tensor usespsi^-4 instead.
Metric/scalar jets require second derivatives; K and A require first derivatives.
There is no claim to computed K Hessians. Binary correction regularity is addressed by the explicit interior plateau
derivation below; physical convergence and measured enclosure remain required.

## Minimal acceptance table (declared before physical runs)

| Distinct obligation | Check and acceptance | Status |
|---|---|---|
| Seed conventions | Schwarzschild closed form; mass scaling; rotated/translated Kerr. Scaled field error<1e-11 at moderate parameters | pass, including production near-horizon checks |
| Derivatives and vacuum | Independent metric/K finite differences at3 steps; expected convergence before rounding; mass-normalized exterior H/M RMS<1e-7 | pass at the sampled exterior, axis, horizon and interior points |
| Boost and causal slice | Generic spin/boost andGamma10; positive metric and slice margin at every evaluated point; E/P/J extrapolation error<1e-4 with two angular levels | sampled checks and charges pass; nonspinning global slicing proved |
| Linearization | One curved manufactured operator case and centered JVP difference sequence; relative error<1e-7 | pass; small device control also passes |
| Moderate binary | Three increasing resolutions; independent fixed exterior near/bulk H/M RMS<1e-6 and max<1e-4 with decreasing errors; positivepsi; charge changes<1e-3 | two sequence grids meet bounds; bulk convergence fails |
| Consumer/horizons | Bound checkpoint and field/gradient roundtrip<1e-12; initial-time expansion RMS<1e-7 at3 angular orders, mass/spin changes<1e-4; modified-region enclosure | isolated spin.99 and finest moderate binary pass |
| QI compatibility | One unchanged saved generic QI fixture, default-path bitwise comparison on identical build/platform | pass |
| Extreme binaries | Separatechi=.99 andGamma10; same independent physical requirements, measured horizon properties andd/Mirr calibration to50 within1% | pending |

Residuals use mass-normalized physical tensors, separately for near-hole,
bulk and modified regions. Failed gates are retained. A documented verifier
roundoff floor is not permission to raise the physical tolerance silently.
Existing source tests may serve as controls but not as fresh integration proof.

## Backend/run policy

Use existing measurements to select one production backend per stage. At most
one small representative CPU/GPU comparison if selection is unresolved; no
backend/Krylov matrix. Record timing/memory from necessary physical runs.
Local seed development uses one CPU thread. Perlmutter numerical work requires
an allocation; single-GPU work uses a fractional shared interactive allocation.
No extreme binary runs before the minimal prerequisites pass. A failed case
is diagnosed before another resolution or continuation is attempted.

## Integrated family dispatch and current evidence

The explicit family API now selects QI=0 or R0=M trumpet=1 for both active
holes. `HiSpID_Config` retains its original ABI and existing entry points
retain QI behavior. The family is carried by the context, host/device cache
builder, sampler and operator hook. Python accepts `seed_family='trumpet_r0_m'`
on creation or as configuration metadata; it rejects unsupported families
and conflicts with checkpoint metadata. Sampling-only contexts remain CPU.
CUDA seed and cache dispatch have passed the small production control below.

The trumpet extrinsic curvature now contracts the rest-chart spacetime
connection with the boosted normal before transforming spatial indices.
This avoids differentiating the inverse boosted ADM metric. On the same
Gamma10 three-point check the momentum RMS falls 1.337e-9, 8.358e-11,
4.979e-12 across step halvings; K derivative error reaches 7.438e-11 before
rounding. Metric/Hamiltonian finite differences still reach a rounding floor.
The original and new records remain separate; the exact pre-refactor source
is retained as `validation/trumpet/seed-restconnection-executed.hpp`.

The integrated API check gives seed/sampler agreement 2.591e-16 and generic
rotation/translation/mass-scaling errors below 6.890e-16. The exact single-seed
weighted collocation residual is 7.307e-16, and centered JVP differences are
below 2.119e-11. No nonlinear binary solve is implied by these controls.
One generic QI default-path fixture, compiled from the pinned old source and
current source with the same compiler, is bitwise identical for physical
metric, curvature tensor and metric gradients at three points.

Version1 checkpoints keep their QI meaning and serialization. Version2 adds
the required `seed_family` field. The AthenaK reader constructs the matching
native sampler and uses the trumpet horizon radius
`sqrt(m^2-a^2)` rather than the QI half-radius. A standalone build using the
actual AthenaK reader gives zero difference from Python/native physical
fields and gradients for v1 and v2. Both readers reject missing family,
unknown family and a v2 payload mislabeled as v1. The first macOS test launch
failed to resolve the native library's relative install name; the corrected
runner uses the native worktree as its explicit cwd and retains the failure.
This is reader/sampler evidence, not a full AthenaK import or horizon run.

Final local API image SHA256:
`439c82359e044269ecab3aa01d2f43f7e3ccae69b85077d062b644d8b356c71e`.
The final source/image receipt binds the retained local image and test results.
Earlier development records are historical, not current-image qualification.

## Puncture regularity: first analytic obstruction

For the unboosted Schwarzschild member, R=r+m, psi=sqrt(R/r),
K=m/R^2 and conformal A^2=(8/3)m^2 R^2/r^6. Linearizing the scalar equation
at fixed free data and zero vector variation gives the leading operator
`Delta u - 11 u/(4 r^2)`. A spherical harmonic of degree l therefore has
regular homogeneous exponent
`p_l=-1/2+sqrt(l(l+1)+3)`, including p_0=sqrt(3)-1/2, about1.232.
This scalar-sector result already permits noninteger puncture powers.
It is not the indicial spectrum of the full coupled spinning/boosted system.
Companion terms in the chosen scalar split also contribute finite constants
which the correction must accommodate.

Consequently no exponential convergence or C2 correction at the puncture
is assumed. The existing basis enforces regularity along axes away from the
excluded foci; endpoint approximation and coupled indicial behavior still
need assessment before a binary resolution sequence is accepted. The family
metadata is separate from the current basis identifier so future basis
changes cannot silently reinterpret existing unknown arrays.

## Independent seed/operator progress (2026-10-05)

A three-level **coarser verifier** sequence resolves the previous roundoff
ambiguity without changing the seed or tolerance. On the original exterior
points, both physical constraints and metric/K derivative differences decrease
approximately16-fold per halving for static, generic spinning/boosted, and
Gamma10 seeds. Finest Gamma10 H/M RMS are9.299e-9/2.139e-8, below1e-7.
The exact executed driver and probe hash are retained in
`validation/trumpet/seed-coarse-refinement.json`; future runs can use
`check_trumpet_seed.py --step-scale 8` with the compiled probe. These exterior
points do not replace the pending high-spin near-horizon physical checks.

The independent ADM observer uses raw physical tensors and Cartesian finite
differences, not native conformal derivatives. Spin.99 and generic spin/boost
extrapolated charge errors are5.49e-10 and3.19e-9 (scaled by max(1,|expected|)).
Gamma10 energy is already accurate at the initial radii, but momentum has
large finite-radius corrections: the initial256--2048 range fails. Moving
outward and reusing completed surface integrals gives a4096--32768 range with
scaled charge error2.78e-6, angular change4.93e-7 and adjacent radial-fit
change1.80e-5, all below the original1e-4 criterion. This checks positive
momentum in the input velocity direction and the Lorentz-transformed angular
momentum of the generic case. Failed shorter-radius fits remain available.
The first charge script's final bookkeeping call used a nonexistent method;
case results were saved before that error. Their immutable executed source and
subsequent image-hash verification are retained; no successful rerun is invented.

A distinct manufactured case uses four sinusoidal scalar/vector fields on the
generic curved trumpet background. The independent observer builds connections
from metric values, forms contravariant longitudinal tensors and differences
their fluxes. Scalar Laplacian, all three vector components, and conformal
Ricci scalar converge at fourth order; finest scaled errors are at most2.12e-9
against the native operator, below1e-7. Together with the earlier centered JVP
control this completes the local operator obligation, not a binary solve.
`validation/trumpet/seed-evidence-receipt.json` binds the retained results.

Production selection uses the existing frozen measurements: at40x80x16,
CUDA HiSpID/GMRES solve median.396s versus1.70s on OpenMP16, with ready times
5.86s versus7.32s. Select CUDA/GMRES for the solve and existing host geometry
construction; use the CPU sampler for the AthenaK consumer. These are historical
selection data, not trumpet timings or acceptance of the optional GPU setup.
A single small device control covers seed dispatch and residual/JVP cache
consistency, including optional execution-space setup; no physical matrix is
repeated. Production builds use row power3, as in the selected measurements.
The local development library uses its original row convention; it is not
substituted silently for the selected production image.

Perlmutter build allocation59402081 uses one GPU and32 logical CPUs in
shared_interactive, scratch root`/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005`.
The archive contains source fromff5e8ae. Build commands are retained in
`validation/perlmutter_trumpet_build.sh`. CPU sampler and CUDA solver images
are separate stages; no old benchmark or baseline binaries are rebuilt.

For the nonspinning boosted target, spacelikeness can also be established
analytically over the entire punctured domain. The rest Schwarzschild trumpet
has alpha=r/R, beta=m*x/R^2 and gamma^ij=(r/R)^2 delta^ij, R=r+m. For speedv<1,
`q-alpha*sqrt(v_i gamma^ij v_j) >= 1-v*m*r/R^2-v*r^2/R^2
=1-v*r/R >0`. Thus both factors of the slice margin
`q^2-alpha^2 v_i gamma^ij v_j` are positive for everyr>0, includingGamma10.
This establishes the nonspinning target's direct-boost slicing, not a global
claim for arbitrary spinning/boosted combinations.

The native CPU sampler and CUDA libraries compiled successfully and allocation
59402081 released. CUDA image`dca1f9b73cbd4cb43c4dd3dfbf1df93688b5cda4b18e3fe4ead82720a6653c74`
passes the one small device control on A100-80GB in allocation59402262:
seed differences<=4.33e-15, host-cache residual/JVP differences<=1.80e-15,
execution-cache residual/JVP differences<=1.69e-13. It uses one8x12x4 generic
binary state and three seed cases, not nonlinear solutions. Necessary control
times are recorded, but their small-grid/warmup ordering does not justify a
large-grid setup speed ratio. The full AthenaK build fromea925d1e with its
pinned Kokkos6739bc6 and the CPU sampler library completed in that allocation.
The allocation released normally. Executable hash begins49459b582a17; full
consumer build records are committed in AthenaK branch codex/hispid-trumpet-pgen
atfbcc1db7. Full import and horizon runs are still pending.

## Binary correction regularity and first physical run

The first moderate binary uses a declared interior regularization, rather than
assuming that the unmodified trumpet correction is C2. With `far_radius=0`
and `inner_flatten=1`, the existing g=0 plateau makes the **correction** equations
exactly `Delta u=0` and `Delta b + grad(div b)/3=0` near each focus. The source
terms vanish, the operator metric is Euclidean and the operator connection and
its derivatives vanish. Bounded solutions have removable puncture singularities
and smooth Cartesian extensions. The actual physical seed and its r^-1/2
scalar are retained. The vacuum equations are only recovered where g=1;
physical acceptance therefore requires measured horizon enclosure of all g<1
regions. The transition's flattened connection is part of the explicitly
modified interior operator, not a claim to an alternative vacuum metric.

This choice supplies regular endpoint conditions for the coupled system even
with spinning/boosted seeds. At infinity the old bounded modal coefficients
retain the correction's 1/r falloff. The mapped modal basis preserves smooth
Cartesian Fourier modes (higher powers fit within its parity-preserving C2
cap), so no new checkpoint basis is required for this regularized construction.
Actual three-grid convergence is still required: suitability is not accuracy,
and no exponential convergence claim is made across the smooth transition.
The earlier fractional-power obstruction remains applicable to **unmodified**
trumpet punctures, not to the flat homogeneous correction plateau.

`examples/trumpet_configs.py` specifies the first generic unequal-mass moderate
case: the previous moderate masses/spins/boosts, no far attenuation, and
inner radii .15/.5 times each seed's smallest laboratory horizon radius.
These radii are only enclosure candidates until measured. A24x48x8 CUDA/GMRES
pilot is diagnostic, with no three-grid or physical-acceptance claim. It stores
the solved iterate before physical sampling, including failed solves.

The exact spin.99 single-trumpet consumer control uses three angular orders
8/12/16, the existing strict expansion/area/spin checks and independent
CPU-consumer migration. Export and migration now carry the chosen family and
place trial horizon points at the trumpet radius; legacy QI defaults retain
the half-radius. The shared off-grid point selector also uses the family radius.

The plateau condition has a focused implementation check, not just an intended
boundary condition: six points inside the two generic spinning/boosted holes'
g=0 plateaus, with nonzero quadratic scalar/vector jets, produce exactly the
flat scalar/vector residuals (absolute difference zero in the retained local
CPU control). This tests the actual coupled equation routine; the metric and
CTT source outside the correction operator are not replaced by a flat seed.

## Completed first physical runs and diagnosed limitations

The isolated trumpet withseedchi=.99 passes the full AthenaK import and
horizon controls at angular orders8/12/16. Finest measured horizon mass is
1.00000000000409, spinSz=.989999999999997 and expansion RMS4.538e-11;
relative area error is5.80e-11. Zero evolution is verified. Separate-process
CUDA-image/CPU-consumer sampler arrays are bitwise identical, including
metric gradients; the local mirror independently verifies their equality.
Allocation59402564 completed and released. No interior source modification
is used in this exact seed control, so enclosure of binary g-balls is separate.

Independent near-axis/horizon/interior physical checks now pass for both
spin.99 andGamma10. The original spin.99 verifier sequence had clear
fourth-order truncation but failed the tolerance; two additional smaller
steps (only that case) reduce H/M RMS to6.147e-9/4.307e-10. The Gamma10
sequence reaches4.447e-8/4.655e-8 without additional runs. The failed coarse
step records are preserved. An extra step initially could not start because
the main job held its reserved communication ports; the independent CPU
observer subsequently ran with MPI disabled and zero reserved ports.

The generic moderate24x48x8 pilot converged in10 Newton/44 Krylov iterations,
with CUDA solve time.162s. It is physically **unaccepted**: near H/M RMS
5.83e-4/3.45e-4, bulk8.95e-6/1.99e-4. Halving the independent observer's step
changes H by at most9.99e-9 and the momentum norm by7.59e-10, excluding verifier
truncation as the dominant error. Nonzero high azimuthal modes remain in the
retained correction. Both physical sampling regions have wholly unmodified
stencils, positivepsi and positive metrics.

The warm-started40x80x12 run also converges internally (9 Newton/41 Krylov,
.429s solve) and remains physically **unaccepted**. Near H/M RMS improves to
1.14e-4/9.43e-5, bulk momentum improves to8.69e-5, but bulk Hamiltonian grows to
1.49e-5. This is not a convergent three-grid sequence. Extrapolated charges are
already stable to4.71e-5 between these two grids, which does not override the
constraint failure. Setup for the refined run uses the selected16-thread
OpenMP host construction and the same CUDA solver; no backend matrix is run.
A directional replay of the retained 40-grid polynomial identifies radial
and polar underresolution as dominant at that resolution; azimuthal replay
errors are much smaller. This diagnostic does not establish the asymptotic
resolution requirements of subsequent grids.

`first-physics/`, `refinement-physics/` and `first-physics-receipt.json` retain
original logs, checkpoints, coefficients, observers and source/image bindings.
The refinement allocation59403087 completed and released. The target factory
now specifies both requested separate equal-mass cases and accepts a measured
component Mirr for the Gamma10 separation update; no target binary has run.


## Further refinement and production horizon backend

The selected CUDA solver and 16-thread host geometry produced these additional
moderate-binary results. All remain diagnostic, with positive sampled psi and
unmodified exterior stencils. RMS tolerances remain 1e-6 for both constraints.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | Solve seconds |
|---|---:|---:|---:|---:|---:|
| 80 x 160 x 16 | 2.597e-5 | 9.208e-6 | 1.048e-7 | 1.267e-5 | 2.345 |
| 160 x 256 x 16 | 2.248e-6 | 1.278e-6 | 2.184e-8 | 2.860e-6 | 16.856 |
| 256 x 256 x 16 | 7.857e-7 | 1.312e-6 | 3.291e-7 | 6.880e-6 | 83.469 |

The attempted 160 x 320 x 20 input exceeded the native per-dimension cap of256
and failed before solving; its receipt is preserved. The runner now checks
this limit and accepts an explicit polar dimension. No failed physical test
has been reclassified as accepted. Refining the radial dimension alone from160
to256 worsens bulk errors, so further blind refinement is stopped. Halving the
independent observer step at256 changes H by at most1.71e-8 and M norm by
1.42e-9; the momentum failure is not explained by observer truncation.
The same-grid replay gives native physical-equivalent residual maxima below
4.3e-13 in the exterior bins. Doubling only the azimuthal sampling of the fixed
polynomial gives near component RMS below9.8e-10 and bulk below5.4e-11. Thus
azimuthal aliasing is much smaller than the off-grid physical errors; further
Fourier refinement or tighter Newton tolerances is not justified. The evidence
points to unresolved radial/polar structure. This is a diagnostic, not a new
nonlinear solve or a physical acceptance test. The first doubled-azimuth
context was rejected by the conservative full-solve memory bound; the replay
then used minimal Krylov controls because it performs no solve.
The remote retained checkpoints are bound by hashes in each result record;
compact physical arrays and result JSON are mirrored locally.

The AthenaK direct native callback now has an opt-in parallel host path, with
a serial cache warmup and worker-exception propagation. The existing spin.99
lmax8 control took12.286s with16threads versus102.03s previously; invariant
mass/spin/area agreement is2.634e-15 scaled and expansion RMS4.420e-11.
This is one matched backend control, not a new angular or performance matrix.
A negative callback test correctly aborts with no horizon data. Its original
verifier incorrectly required the header-only summary file to be absent;
`parallel-consumer/assessment.json` corrects that assertion from retained logs,
without rerunning the numerical case. This changes no physical tolerance.
The selected horizon backend is now this OpenMP consumer. Neither target
binary is yet physically accepted; the moderate convergence gate is unresolved.

## Focused-map trial

The next change uses the existing invertible map family with radial stretch
0.05 and polar stretch 3, concentrating points toward the two punctures. The
physical free data, interior modification and tolerances are unchanged. These
map parameters already have independent Cartesian Hessian and flat-vector
inverse evidence in `private_operator_controls_focus05_k3.json`; the new image
gets one small CUDA/host residual/JVP comparison, not repeated seed tests.
The map-independent trumpet geometry and legacy-QI evidence remain applicable.

CMake now exposes `HISPID_RADIAL_STRETCH` and `HISPID_ANGULAR_STRETCH` and
propagates their definitions to derivative-header consumers. Defaults remain
0.2 and 2. The existing versioned map identifier distinguishes checkpoints.
`run_trumpet_pilot.py --initial-source-library <original image>` explicitly
permits a changed-map initial guess only after checking that image's hash,
the source map identifier, seed family and physical configuration. It reuses
the existing modal remapper and always solves afresh. It does not load source
and target native libraries into one process or inherit physical acceptance.

`validation/trumpet/focused-map/plan.json` and `run_perlmutter.sh` bind the
single 160 x 256 x 16 trial, its source checkpoint and build/run recipe.
The new image reuses the existing pinned CUDA Kokkos package. No binary pass
or adoption of the focused map is claimed before the physical observer runs.


The shared pilot runner also accepts `--case spin99` and `--case gamma10`,
using the already specified separate target factories. For the boost case,
`--measured-mirr <component value>` sets the next separation to50 times that
value; it does not certify the new trial's measured separation ratio. The
moderate default and all physical acceptance requirements are retained.


The focused-map build and selected host-cache CUDA control pass (residual/JVP
relative differences 6.49e-15/1.81e-16). Its 160 x 256 x 16 solve converges in
7 Newton/64 Krylov iterations, with 17.46s solve time. Physical accuracy worsens:
near H/M RMS 3.42e-6/7.40e-5 and bulk 1.21e-6/1.91e-4. **The focused map is
not adopted.** All results and the image hash are retained. The next observer
samples physical tensors at actual collocation nodes, to distinguish
interpolation error from any solver/sampler inconsistency before another
representation change. The preliminary serialization failure is recorded;
the corrected observer retains raw arrays before assembling its JSON report.


At twelve actual collocation nodes nearest the largest retained off-grid
momentum errors, the independent physical observer finds near H/M RMS
4.45e-9/4.87e-11 and bulk 2.75e-9/1.17e-12. All stencils are exterior.
This supports genuine meridional interpolation error rather than a mismatch
between the native equations and sampled physical tensors at those nodes.
The original map is retained. The next single trial increases only the polar
extent from256 to384, holding radial256 and azimuthal16 fixed.

The native/consumer extent limits are now declared together in `HiSpID.h`:
256 radial,512 polar,256 azimuthal. Only the former hard-coded polar cap changes;
the conservative aggregate and device-memory checks remain in force. The
largest permitted integer grid product remains within32-bit index range, and
memory checks reject impractical combinations before constructing workspaces.
The source audit found dynamically sized polar matrices/loops, with no fixed
256-entry buffers. A small actual reader/sampler check at polar384 and an
oversized-polar rejection are scheduled before the full trial. Full AthenaK
horizon qualification of a refined solved checkpoint remains separate.


The polar384 actual-reader control passed: native and consumer tensors and
metric gradients agree exactly, and polar513 is rejected. The refined production binary solve completed under allocation59404717,
which has released. It is the first grid satisfying the near/bulk physical
bounds: H RMS1.63e-7/4.09e-8 and M RMS4.42e-7/6.06e-7; the largest exterior
H/M point residual is1.70e-6. The minimum sampled psi is1.02899686 and all
exterior stencils are unmodified. It used6 Newton/44 Krylov iterations,
126.0s solving and350.4s total. This single grid does not establish the
required three-grid convergence or horizon enclosure.


The next convergence sequence is declared in `polar-refinement/convergence-plan.json`:
192x384x16,224x448x16,256x512x16 on the original map, with both meridional
extents increasing and Fourier resolution held at its independently checked
value. A uniform restart64 keeps the finest conservative allocation bound
within64GiB; it changes a solver control, not the equation or physical gate.
The existing 256x384 result remains a directional control, not a selectively
chosen member of the new sequence. No sequence acceptance has been assigned.

### Declared sequence result (allocation 59405258)

The three prescribed solves completed. The independent assessment is retained
in `validation/trumpet/convergence/assessment.json` and **fails** the declared
convergence gate:

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS |
|---|---:|---:|---:|---:|
| 192x384x16 | 1.128e-6 | 3.772e-7 | 5.764e-8 | 4.419e-7 |
| 224x448x16 | 2.081e-7 | 1.654e-7 | 8.708e-8 | 4.438e-7 |
| 256x512x16 | 1.703e-7 | 1.389e-7 | 2.622e-7 | 2.857e-7 |

The coarse near-H error exceeds the 1e-6 RMS bound; both finer grids meet
all exterior pointwise/RMS bounds. Near errors decrease, but bulk H increases
and bulk M is not strictly decreasing. Adjacent scaled charge changes are
8.06e-8 and 2.28e-7. This is not accepted resolution convergence. The finest
solve used 165.51s, 468.29s total and 10,357,364 KiB peak resident memory.
No thresholds have been relaxed and no alternative subset is substituted.

The full OpenMP AthenaK consumer now builds with the larger polar extent;
its image and source hashes are retained under `convergence/consumer`.
The finest checkpoint remains diagnostic. A separate-process producer/CPU
sampler comparison followed by three-order initial-time horizon measurements
is the next independent obligation, not a transfer of physical acceptance.

The finest checkpoint's separate-process sampler migration passed at all2454
witness points: every recorded field and metric-gradient scaled difference is
zero. The consumer is the pure reference CPU library, with no Kokkos runtime
images loaded. The report is retained in `convergence/horizon256/migration.json`;
raw producer/consumer arrays remain in the report's Perlmutter artifact paths.
Allocation59406161 then loaded that diagnostic checkpoint into the full AthenaK
pgen, reporting ADM-to-Z4c roundtrip relative error4.13e-16 and16-thread horizon
geometry. Horizon convergence and enclosure remain pending; successful import
alone does not qualify them.

Halving the finest checkpoint observer's spatial stencil leaves bulk H RMS
at2.619e-7 (original2.622e-7) and bulk M RMS at2.857e-7. The increasing bulk H
error is therefore not explained by the observer step size. This targeted
replay is retained in `convergence/stencil256.json`; it makes no new solve or
acceptance claim. Near H changes modestly from1.703e-7 to1.726e-7.

A retained-polynomial DCT diagnosis (`convergence/spectrum_modes.json`) finds
vector-component final-eighth coefficient envelopes decreasing across the
sequence, chiefly in m=4 modes. Scalar tails instead grow and are dominated
by m=6 at the finest grid. These are modal P coefficients, not physical errors.
The FFT convention is checked against direct Chebyshev sums on selected lines.
The unweighted collocation residual maxima also grow, despite weighted stopping
below1e-12. A single same-grid tighter-stop control is declared in
`convergence/tolerance-plan.json`; no free data, basis or physical gate changes.
It tests a specific possible source of modal contamination before another
representation change. Its result is not yet known.

### Moderate binary horizon qualification

Allocation59406161 completed the checkpoint-bound CPU migration and AthenaK
initial-time horizon schedule l=8,12,16 plus the l=16 quadrature increase from
32 to48 polar points. All rows and the aggregate horizon/enclosure checks
pass. Results and raw finder logs are under
`convergence/horizon256/surfaces`; the binary data remain diagnostic.

| Component | Horizon mass | Irreducible mass | Coordinate-spin magnitude / mass² | Final expansion RMS | Enclosure margin after refinement buffer |
|---|---:|---:|---:|---:|---:|
| +x | 0.6007648746315 | 0.5887780113536 | 0.3895967639598 | 2.571e-10 | 0.2507675204490 |
| -x | 0.4008630456157 | 0.3935773862669 | 0.3726775137246 | 3.149e-9 | 0.1596480351485 |

Maximum mass and dimensionless spin-vector changes are1.63e-11 and1.01e-11;
maximum area changes are3.34e-11 spectrally and9.88e-12 with quadrature.
The spin is the coordinate rotation integral, not an approximate-Killing-vector
measurement. Continuous enclosure applies to the retained harmonic surfaces
with an empirical refinement buffer, not a rigorous exact-surface error bound.
This closes the consumer/horizon obligation for this checkpoint, independently
of the still-failed constraint-resolution gate. The schedule used113.24,106.63,
154.98 and283.23s on16 CPU threads. The allocation released; the queued
same-grid tolerance control started as59406434.

The same-grid tighter-stop control completed its diagnostics but **failed**
the requested nonlinear convergence: after4 Newton/76 Krylov iterations,
the line search stopped above1e-14 (largest weighted residual9.56e-14).
Bulk H RMS changes only2.622e-7 to2.601e-7; bulk M worsens2.857e-7 to2.976e-7.
Near H/M RMS become1.700e-7/1.320e-7. Extra iterations therefore do not resolve
the failed bulk-resolution trend. The failed checkpoint is retained remotely;
result and raw physical observer are in `convergence/moderate256tight`.
The run took127.05s solving and436.16s total; allocation59406434 released.

A small local differentiation diagnostic compares the current D*D matrix
with an analytic barycentric second-derivative formula on Chebyshev T2/T17
at192/256/384/512 nodes. At512 nodes the T2 scaled endpoint error falls from
1.97e-7 to6.70e-8, while T17 endpoint error is essentially unchanged. Interior
errors are smaller. This identifies a candidate matrix-construction improvement,
not the cause of the binary failure. The diagnostic, raw output, compiler and
platform receipt are retained; Apple arm64 long double is not extended precision.
The attempted remote step found an already-released allocation and ran no work.
No production matrix construction or acceptance criterion has changed.

### Opt-in analytic differentiation matrices

`HISPID_ANALYTIC_MATRICES=ON` now constructs Chebyshev-root derivative matrices
using barycentric off-diagonal identities and trigonometric node differences.
It avoids the cubic D*D product; construction is quadratic. Host construction
uses long double, device construction double, with the same mathematical
entries and row-sum diagonals. Both shared libraries receive the same build
option to keep shared inline definitions consistent. The default is OFF,
retaining the existing path. The interpolant, maps, checkpoint parameterization
and physical equations are unchanged.

The independent modal automatic-differentiation control passes locally
(maximum scaled error2.44e-13). The analytic option is experimental: a single
production-path trial, with its small host/device cache controls, is declared
in `validation/trumpet/analytic-matrices/plan.json`. Its independent physical
errors must improve before considering adoption. No geometry-seed tests or
backend benchmark matrix are repeated.

The analytic-matrix CUDA build passed its prerequisite controls on Perlmutter
(allocation59406945). Independent modal AD error is2.44e-13. Small coupled
residual/JVP comparisons against the host reference give1.80e-15/2.03e-16 for
host caches and1.68e-13/4.29e-14 for device caches. The candidate image is
`ee663d0cf6ff7f5cee51825cf7121754aa10cdc160daea4b71d4c6c2b0ad15cd`.
Build, image/source hashes and control evidence are retained in
`validation/trumpet/analytic-matrices`. The binary comparison is running;
these controls alone do not establish an accuracy improvement or adoption.

The analytic-matrix warm-start trial completed, but the initial checkpoint
already satisfies its1e-12 stopping criterion (largest weighted residual
4.14e-13), so it took **zero Newton/Krylov iterations**. This is a compatibility
result, not an independently recomputed nonlinear solution or proof of improved
accuracy. Independent near H/M RMS remain1.703e-7/1.389e-7; bulk remains
2.622e-7/2.857e-7. The matrix option stays experimental and OFF by default.
The trial and physical arrays are retained under `analytic-matrices/moderate256`.

With the seed, operator, checkpoint and moderate-binary horizon checks passed,
and two moderate grids within physical-error bounds, an exploratory equal-mass
aligned-spin chi=.99 case is now underway. It does not wait on further matrix
microdiagnostics, and does not waive the unresolved bulk-convergence gate.
The declared configuration is160x320x32, separation12, seed mass.5 per hole,
using the established host-geometry/CUDA-GMRES image and original matrices.
See `extreme-pilots/plan.json`; allocation59407208 started only after59406945
released. This case remains diagnostic until its independent physical and
horizon checks and the remaining solver convergence obligations are satisfied.
The separate Gamma10 head-on investigation remains required.

### First aligned-spin chi=.99 pilot

The equal-mass seed-chi=.99 case at separation12 and160x320x32 completed
under allocation59407208. It converges internally in9 Newton/132 Krylov
iterations (54.44s solve,267.73s total), with minimum sampled psi1.02894396.
It **fails** the unchanged independent exterior constraints: near H/M RMS
2.773e-4/7.106e-4 and bulk2.451e-2/1.729e-2. It is not validated high-spin
binary data. Complete result and physical arrays are retained under
`validation/trumpet/extreme-pilots/spin99_160`; the image-bound checkpoint and
unknowns remain in the recorded Perlmutter paths. No horizon properties are
claimed for this pilot.

The separate nonspinning head-on Gamma10 pilot starts at separation25, to be
calibrated against measured component Mirr later, with160x320x8 on the same
production backend. Its configuration and limitations are declared in
`extreme-pilots/gamma10-plan.json`; allocation59407326 began after the spin
allocation released. A fixed-polynomial azimuthal replay of the spin pilot is
scheduled afterward, to distinguish Fourier aliasing from meridional error
before another resolution or map is selected. No further spin solve has been
launched and neither physical acceptance threshold nor failed flag is changed.

A geometric resolution diagnostic now binds the input configurations and counts
nodes across each fixed g-transition shell. At the failed spin160x320x32 grid,
only5 radial-axis projected nodes and12 polar-axis projected nodes span the
transition, with187 actual meridional nodes in the shell. The moderate finest
grid has26--31 radial projections and3928--5294 shell nodes. These are geometric
counts, not a physical convergence test, but they identify a concrete possible
source of the high-spin interpolation failure. Candidate map(.03,3.5) raises
the spin counts to14 radial projections,22 polar projections and1115 shell
nodes without changing the physical attenuation radii or grid size. It is not
yet selected or solved; the pending azimuthal replay will guide that choice.
See `extreme-pilots/transition-resolution-bound.json` and its reproducible driver.

### First head-on Gamma10 pilot

The separate160x320x8 nonspinning Gamma10 case at initial separation25 completed
its diagnostics in372.41s (279.12s solve). It **fails** nonlinear convergence:
9 Newton/4903 Krylov iterations end with the final linear solve reaching2400
iterations and true relative residual0.897. The largest weighted nonlinear
residual is7.20e-12, above the unchanged1e-12 target. Independent near H/M RMS
are4.321e-3/9.981e-3; bulk H/M RMS are4.560e-7/1.778e-5. The sampled psi minimum
is0.99550 and metric minimum eigenvalue0.98212. Positivity alone is not validation.

The full failed result, physical arrays and log are retained in
`validation/trumpet/extreme-pilots/gamma10_160` and its parent directory.
ADM energy extrapolates to9.4900 in this unaccepted iterate; it must not be
interpreted as a validated physical measurement. No component horizons or
measured-Mirr separation calibration have been established for this case.
The same allocation proceeds to the already-declared spin azimuthal replay.

The spin pilot's doubled-azimuth replay completed at160x320x64 without a new
solve. Native-equivalent residual RMS is at most7.50e-9 near the holes and
3.01e-8 in the bulk, far below the original off-grid errors of order1e-2.
This supports meridional under-resolution, rather than missing azimuthal modes,
as the dominant problem. Raw replay evidence is retained under `extreme-pilots`.

A focused-map spin trial is now declared in `spin-focused-map/plan.json` and
its build/run script. It keeps160x320x32 and all physical inputs fixed, changing
only the radial/polar maps from(.2,2) to(.03,3.5). The failed source coefficients
supply an explicitly image-bound remapped guess; a fresh solve and independent
physical observer are required, with no inherited acceptance. It uses the
original matrix construction and the selected host-geometry/CUDA-GMRES path.
The separate Gamma10 result and its failure remain intact.

### Focused-map spin result and targeted refinement

Allocation59407833 completed the spin99 trial on maps(.03,3.5) at160x320x32.
It converged internally in11 Newton/478 Krylov iterations,94.15s solve and
315.16s total. Independent near H/M RMS are3.146e-6/1.438e-5; bulk are
8.296e-5/3.002e-4. Relative to the original map at the same grid, these are
49--295 times smaller, but **all four RMS still fail** the1e-6 threshold.
Minimum sampled psi is1.0289440. No horizon or physical acceptance is claimed.
Compact results and physical arrays are retained in `spin-focused-map/spin99_160`,
with the comparison in `spin-focused-map/comparison.json`.

This improvement justifies one targeted192x384x32 refinement on the same map
and physical data, using the retained solution as an initial guess. The existing
conservative aggregate allocation formula bounds this at61.032GiB (device27.559GiB),
within the64GiB context limit; larger224x448x32 would exceed that limit and is
not requested. See `refine-plan.json` and `refine_perlmutter.sh`. The Fourier
resolution, production backend and all acceptance thresholds remain unchanged.

The spin192 refinement is running as allocation59408144. A single focused-map
Gamma10 trial is queued as59408229 with an afterany dependency, so it requests
no concurrent GPU. Its original160x320x8 grid, separation25, physical free data
and tolerances are unchanged. The original map resolves its transition with
only3 radial-axis projections/67 meridional nodes; the focused map gives8/472.
This geometric evidence motivates the trial but is not proof of physical
accuracy or improved linear conditioning. The already-passed changed-map
operator control is reused. The failed original checkpoint supplies only an
explicitly image-bound remapped initial guess. See `gamma-focused-map/plan.json`,
`transition-resolution.json`, and the run script. Separation calibration using
measured component Mirr and binary horizon validation remain required.

The focused-map CPU sampler and OpenMP AthenaK consumer built successfully on
CPU allocation59408305, which then released. The isolated outputs are
`build-focused-sampler/libHiSpID.so` (SHA65aaa4c8...) and
`build-athenak-focused/src/athena` (SHA6c2e2301...). Both use the same(.03,3.5)
map as the new CUDA producer; original-map consumers are preserved. Reproducible
commands and complete build logs/image hashes are in
`spin-focused-map/build_consumer_perlmutter.sh` and `spin-focused-map/consumer`.
Build success does not establish checkpoint migration or horizon correctness;
those checks remain to be run on the selected target checkpoint. Spin192 is
still solving under59408144; Gamma10 map trial59408229 waits on its release.

### Focused spin checkpoint import and initial horizon attempt

The spin160 focused-map checkpoint passed producer-to-CPU sampler migration
on2454 points with exactly zero differences in every tested field and metric
derivative. AthenaK imports it with ADM/Z4c relative error4.108e-16. Migration
proof and the failed finder attempt are retained in `spin-focused-map/horizon160`.
These checks establish sampling/import compatibility, not physical acceptance.

CPU allocation59408476 failed the firstl8/n16 horizon row after38.26s. Its
initial radius is positive; the iteration trace develops alternating, growing
shape/expansion oscillations before iteration18 requests a nonpositive radius.
The callback label says "initial geometry", but is used on every flow iteration;
this is not evidence that the initial seed radius is invalid or no horizon exists.
A single targeted retry uses flow_alpha=.2 (the update coefficient scales
linearly with alpha), keeping the checkpoint, initial shape, angular schedule,
600-iteration limit and strict1e-7 expansion threshold unchanged. It reuses the
passed migration proof instead of repeating that comparison. Retry allocation
59408628 uses the CPU/OpenMP16 consumer, with script
`run_horizons160_damped_perlmutter.sh`. Results remain diagnostic until obtained.

The first queued Gamma map allocation59408229 was cancelled by its600s pending
wait limit before any compute step ran (Slurm records elapsed0). Replacement
59408624 has a1800s wait limit and the same afterany dependency on still-live
spin refinement59408144. There is no duplicate Gamma solve or concurrent GPU.

### Focused spin refinement: linear failure and radial-tail diagnosis

The192x384x32 refinement completed under59408144 in786.41s total (482.63s
solve), with peak host RSS8497656KiB. It **fails** nonlinear convergence: after
10 Newton/2809 Krylov iterations the last linear solve hits2400 iterations
with relative residual0.506, leaving maximum weighted residual1.402e-12.
Near H/M RMS are2.987e-6/4.657e-6; bulk are2.307e-4/2.393e-4. In particular,
bulk H worsens versus160x320x32. The failed checkpoint/result/physical observer
are preserved, and no larger nonlinear solve is justified by this result alone.

Retained-polynomial spectra show large radial tails dominated bym4 (including
the sine partner index20). Polar final-eighth vector tails shrink from maxima
35.14 to4.935, whereas radial maxima remain490.0 and459.8. These are modal-P
coefficients, not physical constraints; they implicate radial resolution or
conditioning but do not distinguish those causes. `spectrum.json` binds the
checkpoints and diagnostic source. A fixed-polynomial replay separately changes
radial192->224 and polar384->448 with no nonlinear solve, queued as59408748
after Gamma10 allocation59408624. The replay driver now accepts explicit axis
counts to stay inside existing grid/memory caps instead of requiring doubling.

The damped horizon attempt59408628 keeps positive radii and reaches stable
area7.1731520901, butl8/n16 expansion RMS4.4604e-7 still exceeds1e-7 at600
iterations. This failed result is retained under `horizon160-damped`. No area,
mass or spin from this attempt is accepted. The next focused angular sequence
isl12/16/20 with alpha=.2 and unchanged strict tolerance, using the same bound
checkpoint and previously-passed migration proof. It tests angular truncation
without changing initial data or launching a flow-parameter sweep.

### Gamma10 focused-map failure and symmetry diagnosis

Allocation59408624 completed the focused-map Gamma10 trial in287.18s total
(196.20s solve). It **fails** after4 Newton/3500 Krylov iterations; the final
linear solve uses2400 iterations with relative residual0.1835. Largest weighted
nonlinear residual is1.177e-5. Near H/M RMS are3.640e-3/5.242e-2; bulk are
2.619e-5/3.746e-4. The map does not improve this case and is not adopted as
validated boost data. Retained results are under `gamma-focused-map/gamma10_160`.
The unaccepted iterate has spurious transverse angular momentum; its charges
are diagnostics only.

For exactly coaxial nonspinning data, rotational symmetry about the local
puncture x axis implies scalar andb_x have onlym0, whileb_y/b_z have matching
cos(phi)/sin(phi) coefficients, with no azimuthal circulation. Inspection of
the retained modal-P vectors finds forbidden-sector relative maxima0.053--0.305
in the original-map failed iterate, and0.993--1.000 in the focused-map failed
iterate. These auxiliary-coefficient norms are not physical-field error norms.
They motivate an opt-in symmetry-preserving Krylov search, with full unprojected
Newton stopping/line search and independent physical acceptance retained.
`inspect_trumpet_axisymmetry.py`, `gamma-focused-map/axisymmetry.json` and
`spectrum.json` bind this evidence; no coefficients are changed by those scripts.

The spin fixed-polynomial directional replay59408748 completed without solving.
Refining only radial192->224 gives bulk native-equivalent component RMS
(3.759e-5,1.110e-4,3.626e-5,3.395e-5); refining only polar384->448 gives
(1.590e-5,1.001e-5,3.280e-5,2.942e-5). Both directions have unresolved error,
with the radial Hamiltonian andx-momentum contributions larger. This is not
a physical convergence result. See `spin-focused-map/directional-replay.json`.

### Opt-in axial Krylov sector for the nonspinning head-on case

`HiSpID_set_axisymmetric(context,1)` or Python `solve(axisymmetric=True)` now
restricts the linear search to the exact no-swirl axial sector. It is OFF by
default and rejects nonzero spins or noncoaxial local centers/boosts. The shared
projectors operate on modal unknowns and physical-phi equation rows separately,
with the same algebra in reference and Kokkos kernels. Preconditioner outputs
are projected in unknown space; Jacobian outputs and the linear RHS are
projected in equation space. The initial guess is projected explicitly.

Crucially, nonlinear convergence, diagnostics and line search still evaluate
the full unprojected equations. No nonlinear or physical tolerance changes;
no coefficient clipping or residual acceptance transfer. The full stored basis
and checkpoint format are unchanged. Generic spinning binaries and the default
QI/Bowen--York paths are unaffected by the opt-in setting.

The compact control uses an unequal-mass boosted binary, compares unrestricted
reference and restricted reference/CUDA solves, injects forbidden scalar/vector
modes into the restricted initial guesses, and rejects spinning/transverse-boost
configurations. Local reference and Perlmutter CUDA controls pass; full residual
stopping passes in every control. Evidence and commands are in
`validation/trumpet/axisymmetric/{local-control,control,plan}.json` and
`run_perlmutter.sh`. These are algorithm-equivalence controls, not independent
physical validation of the coarse test binary.

Allocation59409143 passed the CUDA control and started the same160x320x8
focused-map Gamma10 case with this option. It uses the same original-map source
checkpoint, explicitly remapped, with identical physical data and tolerances.
The symmetry option is a candidate; no Gamma10 improvement is claimed yet.

The spin160l12/n24 horizon row passed for both components: mass0.500017851906,
Mirr0.377763585531, coordinate chi0.989929520944, expansion RMS1.694e-8. Retained
surface lower bound0.067529881 exceeds the modified-ball radius0.035266840.
Thel16/20 and quadrature checks remain running under59408772, so these are
preliminary measurements and not angular-convergence or initial-data acceptance.
The completedl12 artifacts and partial metadata are retained under
`spin-focused-map/horizon160-angular`.

The symmetry-restricted Gamma10 trial59409143 completed successfully in the
**full unprojected nonlinear equations**:10 Newton/451 Krylov iterations,
maximum weighted residual5.570e-13,43.55s solve and132.83s total. The prior
unrestricted same-map attempt failed after3500 Krylov iterations and196.20s
solve. This is a convergence improvement, not a matched-accuracy benchmark.
Angular momentum components return below1e-13 and net momentum below2.4e-12.
However, independent near H/M RMS remain3.640e-3/5.242e-2 and bulk
2.618e-5/3.741e-4, essentially unchanged from the failed same-map iterate.
Thus forbidden-mode contamination explains a solver/charge problem but **not**
the dominant off-grid physical errors; the spatial accuracy problem remains.
The result, arrays and comparison are retained under `axisymmetric/gamma10_160`
and `axisymmetric/comparison.json`. No Gamma10 physical or horizon acceptance
is claimed, and no stopping threshold was changed.

The converged axial Gamma10 checkpoint was re-observed at half the finite-
difference step: near H/M RMS3.63954e-3/5.24243e-2 and bulk2.61808e-5/3.74101e-4
remain essentially unchanged. The observer step is not the dominant error.
Fixed-polynomial replay59409707 separately refines radial160->192 and
polar320->384. Near native-equivalentx-momentum RMS is6.443e-3 versus1.898e-5
(about339 times larger in radial replay); bulk is2.055e-3 versus2.016e-5.
No new nonlinear solve was used for either diagnostic. Evidence is retained
in `axisymmetric/stencil.json` and `axisymmetric/directional-replay.json`.

This motivates a single radial refinement to256x320x8, retaining maps, physical
data, symmetry option, polar/azimuthal resolution and all tolerances. It is
running under59410256. Commands and decision are in
`axisymmetric/refine_radial_perlmutter.sh` and `radial-plan.json`.

The256x320x8 radial refinement59410256 converged in10 Newton/386 Krylov
iterations,161.34s solve and272.65s total. Near H/M RMS improved to
8.226e-4/2.421e-2; bulk to4.464e-6/1.272e-4. The reductions are respectively
4.42/2.17 and5.87/2.94 times, but all four RMS remain above1e-6. This supports
radial under-resolution while leaving substantial spatial error unresolved.
The final weighted residual is7.678e-13; no tolerance or physical parameter
changed. Results and comparison are retained under `axisymmetric/gamma10_256x320`
and `radial-comparison.json`.

The next targeted Gamma10 experiment59410842 changes only the radial map
from.03 to.003 at256x320x8. The same fixed interior transition gains35
projected radial nodes versus13, while polar nodes remain18. This is a
resolution estimate, not accuracy evidence. The isolated image must pass a
small axial equivalence control before the solve. Plan and commands are in
`axisymmetric/refocus-plan.json` and `refocus_radial_perlmutter.sh`.

Spin horizon cold-start orders12 and16 both pass, with mass.500017851906
and coordinate chi.98992952095; their retained surfaces enclose the modified
regions. The cold-start allocation59408772 was intentionally stopped after
saving both orders, once new shape-reuse allocation59410814 reproduced them.
The latter uses retained harmonic coefficients solely as initial guesses,
with fresh geometry and every original acceptance check. Times fall from
670.69/1226.91s to35.49/50.42s at orders12/16, with mass/spin agreement within
1e-12. Higher-order and quadrature checks remain pending. This is horizon
measurement on diagnostic data, not physical acceptance of the spin binary.

The spin160 horizon schedule59410814 is now complete and passed: spectral
orders12/16/20 and independent quadrature20x60, with strict expansion criteria
at every order. Final component mass=.50001785189996, Mirr=.37776358547565,
coordinate chi=.98992952098846, expansion RMS=5.976e-9. Maximum mass-relative
and spin-vector changes are1.106e-11 and3.788e-11. The final continuous retained
surface enclosure margin after empirical refinement buffer is.03226304058.
All retained surface, input-guess and log hashes were checked after local copy.
Evidence: `spin-focused-map/horizon160-warm/surfaces/binary.json` and adjacent
artifacts. The allocation exited normally. These successful measurements do
not change the failing independent physical constraints of the input data.

The Gamma10 radial-map control passed (image d51f846821ea7cdd5cd14399090634072ff4d157d0a16e4647bc9aaad3e1a26f);
the targeted solve remains running under59410842. Largest fixed-observer errors
from the preceding Gamma10 and spin checkpoints, explicitly separated from
modified-region errors, are retained in `axisymmetric/error-localization.json`.

Gamma10 stronger radial focusing59410842 completed, but is rejected. At the
same256x320x8 grid, map.003 converged in11 Newton/441 Krylov iterations,
187.33s solve and308.10s total. Independent near H/M RMS are.20795/1.39340;
bulk6.362e-4/5.585e-3. Both greatly worsen relative to.03. Radial final-eighth
modal tails rise from2.35e-7 to1.54e-5 (scalar),2.30e-5 to2.02e-4 (axial
vector), and4.25e-5 to1.72e-4 (transverse vector). Node counts in one transition
are therefore insufficient to select a map. All failures are preserved under
`gamma-radial-focused`; no new physical acceptance is claimed.

The local spectrum witness encountered platform matrix-multiplication warnings.
Its independent direct sum now uses explicit einsum summation and rejects
nonfinite inputs, transforms or witnesses; the original and corrected records
are retained separately. Corrected FFT/direct checks pass for both actual
checkpoints, without changing their polynomial spectra.

For Gamma10 horizon measurement, the existing isolated-seed history records
slow/nonconvergent high-order spherical-harmonic searches of the Lorentz-flattened
surface. A prospective improvement is a per-hole affine spatial chart, stretching
the boost direction for the finder while pulling back gamma, K and spatial
metric derivatives consistently. This changes no slice or initial data. It
requires a separate tensor-transformation/invariant-area control, chart-aware
seed guesses and enclosure, and explicit spin-coordinate conventions before
application to binary data. It is not yet implemented or validated.

## Affine Gamma10 horizon validation

The optional affine spatial chart x=c+J(y-c), with J contracting by1/Gamma
along the seed velocity, preserves the laboratory slice. The consumer pulls
back gamma/K/dgamma and uses J^{-1} E_i J for laboratory coordinate-spin
integrals. Default finder behavior is unchanged. Standalone derivative,
area-element and spin-integrand invariance checks passed locally and on
Perlmutter; full consumer build59411280 passed.

The exact isolated Gamma10 control59411566 passed at orders4/8/12. Max relative
area error7.22e-11, shape error1.13e-10, expansion RMS4.78e-11 and horizon-mass
error3.61e-11, with23.25/59.52/115.11s runs. Producer/CPU sampler migration is
exact. Earlier input-declaration, iteration-cap and report-serialization failures
are retained under `affine-horizon/exact`; only `horizons-v4` passes.

The affine exact sphere is representable at every tested order. Qualification
therefore requires every order to pass the unchanged physical/shape bounds plus
stable mass/spin, rather than monotonic decrease of a nonlinear stopping residual.
The laboratory-chart rule is unchanged. This is not a relaxation of area, shape,
expansion, mass-stabilization or binary constraint tolerances.

The saved Gamma10 binary256x320x8 (map.03) is now being measured using the same
finder, a checkpoint-bound sampler proof, three orders plus quadrature, and
warm guesses reused only after passing rows. No physical binary acceptance is
inherited. A separate feasibility note records that radial512x320x8 is estimated
to fit the existing64GiB aggregate guard; it is not launched or validated.

The first binary Gamma10 affine search59411746 imported exactly (ADM/Z4c
roundtrip4.58e-16) but failed at l8: expansion RMS1.8313e-5 remained above1e-7
at1000 iterations despite stable area. No horizon mass from that attempt is
qualified. The same checkpoint, executable and sampler proof are being used
for l16/24/32 plus quadrature under59411890; no tolerance changed.

The radial cap is now512 in the native header and diagnostic pilot launcher,
with all allocation guards unchanged. A512x16x8 manufactured mapped-operator
control on m0/m1 scalar/axial-vector fields passed: max normalized conformal
operator error5.486e-10 against the existing1e-8 limit. The independent Cartesian
Hessian oracle covers cosine/sine partners. This isolates the changed radial
extent without repeating the full existing mode/backend suite. The test driver
now supports independent polar extent and selected modes/components.

The production512x320x8 Gamma10 solve (map.03, axial GMRES, restart64) is running
under59412030, with a fresh output `radial512/gamma10_512x320-v3`. The analytic
allocation estimate is63.74GiB aggregate and28.10GiB device within the existing
64GiB budget; runtime guards still apply. Two prior launch errors (Python cap
and existing empty output directory) are retained, and neither performed the
large solve. Existing consumer images still have the old radial cap; loading
future512-point data requires a matching sampler/header rebuild. No higher-
resolution physical accuracy or consumer acceptance is claimed yet.

The binary256 affine l16 search59411890 also failed: after1000 iterations,
expansion RMS1.8834243e-5, essentially the same floor as l8. Its area was stable,
but that does not qualify a horizon. The script correctly stopped before
orders24/32 and quadrature. Logs and partial receipt are retained in
`affine-horizon/binary256-angular`. Further angular-only searches on these data
are deferred until the finer elliptic result can distinguish data resolution
from finder truncation; no mass/separation calibration uses these failed surfaces.

The matching radial512 CPU sampler and OpenMP AthenaK consumer built successfully
under59412091. Sampler SHA4c602e408731923c9f426033186e03150c12875fe56ed39bda3bd4afd7e4703e;
consumer SHAde410d82ce2d935cac939cfbc05a9e82747508b9b427ab09e8586262f9b76bf8.
Sources and build directories are isolated from earlier images. A checkpoint-bound
producer/sampler comparison is queued under59412182 after the successful end of
solve59412030; it transfers sampling compatibility only, never physical acceptance.

### Moderate bulk Hamiltonian localization

A targeted CPU-sampler replay on compute step59412030.2 compared the retained
256x512x16 moderate checkpoint with the same polynomial after zeroing only the
scalar cosine/sine m6 rows (6 and14). All physical free data and the18 existing
bulk observer points/stencil sizes were held fixed. The unchanged replay matches
the original H and M diagnostics exactly, using the already-qualified polar CPU
sampler10b346af… and checkpointc1d55c98…. The saved source and result are
`validation/diagnose_trumpet_scalar_mode6.py` and
`validation/trumpet/convergence/scalar-mode6-probe.json`.

Removing scalar m6 reduces bulk H RMS2.6219737e-7 to5.1987097e-9 (50.44x), and
H maximum1.0709601e-6 to1.8493774e-8. Momentum RMS stays2.8570e-7. Thus this
specific retained mode accounts for most of the finest-grid bulk Hamiltonian
error, consistent with the earlier growing scalar-m6 coefficient tail. This
is a sensitivity diagnosis only: the modified polynomial is not a new solve,
not checked in all regions, not exported as initial data, and not accepted.
It does not prove the true solution has zero m6. The next mathematical target
is the mode's representation/conditioning and axis regularity; further uniform
refinement or simply discarding m6 is not justified by this result.

The radial512 production attempt59412030 reached its30-minute allocation limit
without returning from the nonlinear solve. Its preliminary result records
73.11s setup and successful initial-guess transfer, but no final Newton/Krylov
history, solved array, checkpoint or physical diagnostics. Slurm reports TIMEOUT;
the partial result and accounting are retained under `radial512`. This is a
runtime failure, not evidence of nonlinear or physical convergence. The dependent
CPU migration59412182 was cancelled after Slurm marked DependencyNeverSatisfied.

An opt-in-axisymmetry optimization now retains only the m0/m1 preconditioner
factor families, four groups instead of ten at nphi8. Higher temporary output
modes alias existing factors and are discarded by the existing projection;
active equations, full nonlinear residual acceptance and default configurations
are unchanged. The small local contaminated-guess/unrestricted-solution control
passes with the same4.73e-11 maximum scaled field/gradient difference as before.
The CUDA build/control is running under59413039; no production speedup is yet
claimed. Memory guards remain conservative and unchanged.

CUDA allocation59413039 completed successfully. The reduced-factor image
39904a219550a69d54461fb88456b07f0b3c138801f2343008f82363e255dfc6 passes the
same small CPU/reference and CUDA axial equivalence control, including excluded
initial modes and rejection of unsupported spin/transverse boosts. Complete
build/source hashes and control records are retained in `axial-factors`. This
qualifies the factor reduction for the next production attempt, not an accuracy
or speedup claim for the unresolved radial512 physical case.

The next radial512 Gamma10 production attempt is running as59413586 with the
qualified reduced-factor image, identical physical configuration/maps and the
same retained radial256 initial checkpoint. It uses a fresh output directory
`axial-factors/gamma10_512x320-v1` and a one-hour single80GBGPU shared allocation.
The longer limit addresses the observed runtime failure; no convergence bound
or memory guard changed. Consumer migration will be scheduled only after an
actual solved checkpoint is available. The cancelled timeout-dependent job is
not reused.

### Axis regularity diagnosis and isolated C4 experiment

For a smooth Cartesian scalar, mode m6 starts at rho^6 away from punctures.
The retained C2 basis permits a rho^4 m6 term. With t=a^2 and half-separation b,
its coefficient in the orthonormal real Fourier basis is
`-(1-t)^5 P6/(8 b^4)`. Smoothness therefore requires P6 to vanish on t=0 and
eta=+/-1, away from the punctures; those boundary zeros are not imposed by the
capped basis. Barycentric evaluation of the saved polynomial on nine fixed axis
locations gives maximum cosine/sine coefficient magnitudes4.89e-8,5.19e-7 and
8.92e-7 across the three moderate grids. This is an18.2x increase, not convergence
to smoothness. Raw results/source are `convergence/axis-regularity.json` and
`validation/diagnose_trumpet_axis_regularity.py`. The first probe attempt only
failed a metadata-key lookup before writing output; the corrected probe verifies
the retained vector hash and reads the basis from the checkpoint header.

`HISPID_REGULARITY_CAP=6` is now an isolated experimental build option. It uses
r=m through m6, then odd/even caps5/6, guaranteeing C4 rather than only C2 on
axes. Default cap4 and all old checkpoint identifiers remain unchanged. The
new `modal_P_C4prolate_map_v4_r..._k...` identifier distinguishes the continuous
basis while retaining checkpoint text format2. Native AthenaK readers compare
the full identifier and therefore reject C4 files when linked to C2 libraries.
Python metadata and same-basis modal prolongation understand C4; the pilot
explicitly rejects cross-regularity-family initial guesses, even with an explicit
source library. No conversion or acceptance is inferred from matching maps.

The local independent Cartesian manufactured-operator check covers m5/6/8,
cosine/sine partners, scalar and axial-vector components. At24 nodes m5/6 pass,
but m8 has3.29e-8 normalized error and fails1e-8. At40 all pass, max8.78e-13.
Both records are retained in `c4-experiment`. A small exact checkpoint roundtrip,
constant-m6 prolongation and invalid-version rejection also pass. These are
prerequisites only. The CUDA image is building on compute step59413586.3;
its production-path control and one moderate192x384x16 trial must complete
before judging whether stronger regularity improves the binary physics.
The active Gamma10 run retains its already-qualified C2 image and is unaffected.

The isolated C4 CUDA build completed on compute step59413586.3. Its source
manifest matches the working tree. Job59414246 is queued after the active boost
job, so only one GPU is used: it first runs the independent m5/6/8 operator and
metadata control on CUDA, then (only if that passes) a zero-start moderate
192x384x16 solve with unchanged physical inputs and acceptance bounds. No C2
checkpoint is reused as C4 coefficients. Default production remains C2 pending
this actual binary comparison.

### Batched modal transfer construction

A live snapshot during59413586 found the allocated GPU idle with4647MiB resident;
this is one phase observation, not a whole-run utilization profile. Inspection
identified a costly host operation: each radial-block transfer A^{-1}diag(U)
was built by n separate strided vector LU solves. The opt-in
`HISPID_BATCHED_MODAL_TRANSFER=ON` instead forms P diag(U), then uses lower/unit
and upper/nonunit matrix-RHS triangular solves. The Schur recurrence, LU pivots,
preconditioner matrix, PDE and nonlinear stopping are unchanged. Default is OFF.
The shared helper retains the original path and avoids repeating factor logic.

An independent pivoted dense-matrix control at64 and512 rows compares both
transfers and checks A*T=diag(U) against the original matrix, with explicit
finite checks. At512, transfer disagreement is5.11e-16 and independent scaled
residual is7.41e-16 on Perlmutter. The necessary target-node kernel comparison
was0.555839s column-wise versus0.107008s batched (5.19x); local timings were
0.128710s versus0.0321041s. These are single-matrix diagnostics on allocated CPU
cores, not a matched-accuracy whole-solver benchmark. No solve-level speedup
is claimed. Test/source/command records are under `batched-transfer`.

The small axial solve control passes locally. Separate processes loading the
old and batched images produce field/metric-gradient witnesses differing by
at most1.303e-16 scaled, against1e-12. The existing control can now retain its
last axial field witness for direct cross-image comparison; it does not rerun
another physical suite. The matrix test is registered in CMake and Makefile.
Neither the active boost image nor queued C4 image was changed. CUDA compilation
and adoption of this option remain future work if those runs need continuation.

A short read-only debugger snapshot of the running59413586 solver confirms an
active CPU thread inside `cblas_dtrsv -> gsl_linalg_LU_svx -> ModalBlock::factor`;
the debugger detached cleanly. The stack is retained in
`axial-factors/stack-snapshot.txt`. This identifies the transfer-construction
phase directly, rather than inferring it solely from GPU utilization.

The batched C2 CUDA image built successfully. Its first control failed before
solving because the isolated source snapshot had an older Python example module.
Only the queued experimental workflows were refreshed with the current Python,
example and validation sources; the active production workflow was untouched.
The failed import log is retained. The corrected control passes with image
fdaf287805d5f66dc66a72cd9b445188bfa1a1834d818c4eacfc270fedbd1b61;
records are in `batched-transfer/cuda/control-v2.json`. The queued C4 pilot also
passes an import-only preflight, and its refreshed workflow source manifest is
retained as `c4-experiment/workflow-manifest-v2.sha256`. This refresh changes
neither native image nor mathematical input configuration.

The checkpoint-bound migration and affine-horizon command for the active boost
attempt is prepared in `axial-factors/consumer_perlmutter.sh`. It uses the rebuilt
radial512 CPU sampler/consumer, new output directory and the actual producer
image; it is not launched before a solved checkpoint and reviewed diagnostics
exist. No physical accuracy or whole-solver timing conclusion is added here.

The qualified batched image also has a prepared radial512 command in
`batched-transfer/solve_perlmutter.sh`. It retains the current physical inputs,
grid, tolerances and initial radial256 checkpoint, pins the tested native image
hash and writes a fresh output directory. This is a contingency for a terminal
runtime failure of59413586, not a second concurrent run or a performance sweep.
At the latest observation59413586 remained RUNNING at44:23 with no solved
checkpoint;59414246 remained queued behind it. Neither job was interrupted.

### Completed radial512 boost refinement

Job59413586 completed in55:10 and produced checkpoint
`d90bc47d886c395a3acf14e4beb4cd3f5cc72882e948e2682a9e073900b4933b`.
The nonlinear solve converged in10 Newton/352 Krylov iterations,3119.562s;
total workflow3294.567s and peak host RSS5856828KiB (5.59GiB).
Independent near H/M RMS are1.01625e-4/6.12380e-3; bulk1.06906e-6/3.75054e-5.
Relative to radial256 these improve by8.09/3.95 and4.18/3.39, respectively,
but still fail the existing physical bounds. Minimum psi0.99550695, minimum
metric eigenvalue0.98214848, and exterior stencils remain unmodified.
Extrapolated energy9.49075644 changes relatively5.48e-6 from radial256;
net momentum/angular momentum remain at about1e-12 or below.
Results, physical samples and comparison are retained in `axial-factors`.

The batched duplicate is unnecessary and was not submitted. The diagnostic
radial512 consumer was submitted as CPU job59415301 using the prepared
checkpoint-bound migration and affine-horizon script. Its outcome is pending;
no horizon acceptance or boost separation calibration is inferred here.
The dependent C4 production-path job59414246 has started on its single GPU.
Its CUDA manufactured/metadata control passed, with largest normalized error
8.77558e-13 and exact checkpoint roundtrip; the moderate192 solve is running.
The result is retained as `c4-experiment/cuda-operator-control.json`.

The C4 moderate192 workflow completed but the solve failed its second Newton
line search (31 total Krylov iterations,18.852s solve,202.761s workflow).
Both linear solves met the requested0.1 forcing (0.07662 and0.07208 relative
residuals); physical bounds also fail. The first accepted iterate is retained
in checkpoint0f6ae4f8880fc7778a52d056c6cd7346899d1c0c9e8e560a9770bed6ea99aa90.
This does not establish whether C4 improves the converged solution.

One targeted inexact-Newton diagnostic resumes that same C4 checkpoint with
linear forcing0.001, otherwise identical inputs and native image, as single-GPU
job59415408. The pilot now exposes `--linear-rtol` (default0.1 unchanged) and
saves linear history before physical postprocessing, so a later timeout cannot
hide the completed solve's history. The separately named remote driver preserves
the executed original source. CLI import/help and Python syntax checks pass.
No new operator or backend matrix is scheduled.

The radial512 consumer job59415301 ended FAILED after5:05: sampler migration
passed and ADM/Z4c import error was4.39116e-16, but the first l8 affine horizon
reached1000 iterations with expansion RMS9.48912e-6 (required1e-7).
Attempt area12.5655931083 is diagnostic only; no qualified component mass or
enclosure result was produced. The driver correctly stopped before l12/l16,
the second component and separation calibration. Compared with the retained
radial256 l8 expansion1.83132e-5, the floor decreases by about1.93x but remains
unacceptable. Raw surface outputs, migration receipt and failure log are
retained under `axial-factors/consumer-v1`; producer sampling witnesses remain
on Perlmutter under that directory's `raw/` subdirectory.

The C4 tighter-forcing retry59415408 also failed the first resumed Newton line
search:45 Krylov iterations reached true linear relative residual9.07867e-4,
but no trial was accepted. It completed postprocessing in200.466s and is retained
as `c4-experiment/moderate192-linear001`. Tightening linear forcing alone does
not fix this failure.

An opt-in `HISPID_TRACE_NEWTON=1` diagnostic now records the actual linear
relative residual, descent slope, each existing trial's nonlinear reduction,
directional Taylor remainder, minimum psi and nonfinite count. Default behavior
is unchanged; no additional trial is accepted or generated. A small two-step
local control has bitwise-identical iterates with tracing off/on. The retained
C4 iterate is replayed for one step on CUDA; this is not a new physical solve
or a convergence claim. The first two remote builds stopped before numerical
work: omitted CMake test sources, then macOS sidecars in the C-source glob.
The isolated source copy was corrected; both failures are retained in
`c4-experiment/newton-trace`, and job59415685 runs the corrected build/replay.

The existing boost coefficients were inspected on compute step59415685.1
without solving. Radial256-to512 final-eighth modal-P tails change from
(2.35e-7,2.30e-5,4.25e-5,4.25e-5) to
(7.73e-9,2.86e-6,2.27e-6,2.27e-6); radial512 polar tails are
(2.55e-9,1.69e-7,7.41e-7,7.41e-7). These are auxiliary-coefficient tails,
not independent physical errors; both meridional directions remain relevant.
Raw envelopes and transform witnesses are `axial-factors/radial512-spectrum.json`.

The traced C4 replay completed and reproduced the line-search failure. The
linear residual ratio is9.71539e-4 and normalized descent slope-0.999999056;
the modal-P step infinity norm is3.34603e12. All11 tested damping factors keep
psi positive (minimum0.6294 at a full step), and all residual entries are finite.
The nonlinear residual ratio decreases from2.93039e6 at damping1 to1.073274
at the final permitted1/1024 step. The directional Taylor remainder decreases
approximately linearly at small damping, consistent with a large quadratic
nonlinearity. This identifies insufficient tested damping as a concrete next
hypothesis; it does not yet prove that smaller steps converge.

`HISPID_NEWTON_BACKTRACKS` is a build option, default10 preserving the original
trial sequence and cutoff. An isolated24-halving build resumes the same C4
iterate with the existing strict residual and physical bounds; the old source,
images and failed results are untouched. Command is
`c4-experiment/newton-trace/backtrack_perlmutter.sh`. Default local build and
shell syntax checks pass. No new physical acceptance is claimed.

The extended-backtracking job59415954 completed: the first resumed step was
accepted at1/2048, but the second linear solve exhausted2400 iterations at
relative residual0.03459 (target0.001). Totals were2 Newton/2445 Krylov,
245.303s solve and424.004s workflow. Thus deeper damping resolves the immediate
cutoff failure but does not give convergence. Records are `c4-experiment/backtrack-v1`.

Two local preconditioner-only scaling prototypes were screened on the generic
12x24x16 C4 problem. Limiting inverse axis powers to3/4 failed after2689 Krylov;
a1e-12 floor on the dimensionless axis factor failed after3359. Both leave
residual/JVP witnesses bitwise unchanged, but neither reaches the original
nonlinear tolerance. The unmodified C4 baseline converges in6 Newton/135 Krylov.
Both prototypes are rejected and removed from production source; patches,
build options and results are retained in `c4-experiment/scaled-inverse`.
No GPU trial or physical acceptance is assigned to these changes.

The successful coarse baseline motivates continuation rather than another
preconditioner parameter trial. Single-GPU job59416379 uses the existing
extended-backtracking image and unchanged preconditioner to generate a12x24x16
C4 checkpoint, then initializes the same192x384x16 problem from it only if the
coarse nonlinear solve converges. Basis-aware prolongation and all physical
acceptance bounds are unchanged; coarse physical acceptance is not inherited.

The continuation job59416379 completed in7:26 (scheduler exit0), but the fine
nonlinear solve failed. Coarse12x24x16 converged in6 Newton/133 Krylov with
near H/M RMS0.00814/0.00380 and bulk0.00385/0.0220, already outside physical
bounds. Fine192x384x16 accepted its first step at1/1024, then exhausted2400
linear iterations at relative residual0.01810 versus0.001. Total2 Newton/2526
Krylov,252.723s solve,429.353s workflow. Final near H/M RMS0.55864/0.15544
and bulk7.34714/3.50287 are not accepted. Logs and both result records are
retained in `c4-experiment/continuation-v1`; larger artifacts remain in the
matching Perlmutter directory. No job remains running from this experiment.

Source inspection confirms the continuation selects C4 modal-P interpolation
with matching image, basis and map metadata, rather than physical-Fourier
interpolation. This is not yet a numerical verification of transferred fields.
The next diagnostic should compare unchanged coarse and prolonged fields
before any Newton step, and localize high-mode corrections; the final failed
iterate alone cannot distinguish transfer error from nonlinear-step damage.
Do not repeat this coarse-to-fine solve or infer a converged C4 physical result.

A no-solve transfer check now compares coarse and prolonged physical fields
and metric gradients at all102 existing near/bulk/modified observer points.
The same C4 CPU sampler evaluates both representations: maximum scaled field
error is1.64e-13 and gradient error3.19e-12. This rules out a gross transfer
defect at those points, not all possible interpolation errors. The source
checkpoint hash is9738e21d5391e7320b9c4dc40359f912764b6a21f021bbb1de578ae2e7782493.
A separate decomposition of the retained accepted step identifies m6 as the
largest sampled field perturbation: its auxiliary coefficient change reaches
2.083e9, physical correction0.0266, Kij0.505 and metric gradient0.572. Modes
above4 are strongly amplified. This is evidence for investigating C4 high-mode
conditioning/source precision, not a justification to filter accepted data.
No mode-filtered field is saved as a solution and no new physical acceptance
is claimed. Scripts are `diagnose_trumpet_transfer.py` and
`diagnose_trumpet_step_modes.py`; results are `continuation-v1/transfer.json`
and `step-modes.json`. Both use the existing local C4 sampler image with the
same explicit basis/maps; their image hashes are recorded.

Source-precision probe59416928 completed on a Perlmutter CPU allocation,
without an elliptic solve. `diagnose_trumpet_source_precision.cpp` uses the
moderate free data and64 azimuthal samples on cylinders at x=0,2.7,3.5 with
rho=0.1 down to0.0001. It compares the production cached scalar source with
the same geometric jets assembled in64-mantissa-bit long double, recomputing
curvature/laplacian before their existing double casts. At x=2.7,rho=0.001
the m6 Fourier amplitudes are2.68e-16 cached versus1.44e-18 extended; at
x=3.5,rho=0.003 they are6.44e-17 versus2.32e-19. At larger rho the extended
coefficient follows approximately rho^6 before reaching its own roundoff
floor. This demonstrates a scalar-source angular noise floor; it does not
prove that this is the sole cause of the failed coupled Newton update.
The mode-dependent inverse scaling contains `(a*sqrt(1-eta^2))^-r`, so C4
can amplify such nonregular numerical components. Merely assembling the
source in long double cannot be assumed to fix the finest-grid problem.

Records, hashes and reproduction allocation script are in
`c4-experiment/source-precision`. Source changes for a future remedy should
address regular Fourier source evaluation near the axis and stable separation
of background/correction terms, with a derivation and a focused unchanged-
equation control. No source clipping, modal deletion, or acceptance relaxation
has been introduced. Momentum-source precision remains unmeasured by this
scalar-only probe.

An opt-in `HISPID_STABLE_SCALAR_SOURCE=ON` now combines the background
Hamiltonian source in the geometry scalar type before narrowing, then
evaluates correction terms with factored differences of powers. Writing
psi=psi0+u and A=M+L, it uses S0-u*R/8-K^2*(psi^5-psi0^5)/12
+(2 M:L+L:L)/(8 psi^7)+(M:M)*(psi^-7-psi0^-7)/8. The two power differences
are evaluated without background-sized subtraction. The PDE and analytic
JVP are unchanged algebraically; finite-precision constant-source values
change intentionally. Momentum assembly remains unchanged. The option is
OFF by default; its16-byte-per-point cache extension is included in the
Kokkos aggregate memory bound.

Local12x24x16 C4 control converges in6 Newton/129 Krylov (baseline135);
residual witness changes1.56e-16, JVP witness is bitwise unchanged, and a
centered directional difference reaches relative L2 error1.00e-11. Converged
six-point physical fields/gradients change by1.03e-7 in scaled maximum norm;
this is reported as sensitivity, not a strict field-equivalence pass or
physical validation. Records and checker are in `c4-experiment/stable-source`.
A single-GPU shared job59417000 builds the same opt-in path and replays one
Newton step from the original192-grid C4 iterate. It is bounded diagnostic
work, with the existing preconditioner and24-halving option; no full solve
or new physical acceptance is implied.

The stable-source CUDA build59417000 succeeded, but its Python replay stopped
before numerical work because the system interpreter lacks future annotations.
The reproduction script now explicitly loads cray-python/3.12.12; all build
and failure logs are preserved in `stable-source/gpu-v1`. Retry59417104 reused
the image and completed in1:32. Its one Newton step used42 Krylov iterations,
true linear residual9.62332e-4 and14.676s solve time. The step maximum remains
1.63204e12; first accepted damping is1/2048, with residual ratio0.999574880.
Thus scalar-source rearrangement alone does not cure the large-update behavior
and does not justify a new full solve. Results are `stable-source/gpu-v2`;
image316cc298e84f39c0d26fdabce6068f09d088fbd12f71ce7e165f494decb89e87.
It remains opt-in, with no production or physical acceptance. No allocation
from these two jobs remains active. Further work should address representation
and preconditioner conditioning, retaining exact axis regularity rather than
repeating precision-only or damping trials.

A drift-stencil audit finds a separate preconditioner issue: the centered
discretization of a*d2+b*d has negative neighbor weights when |b|h/(2a)>1.
For r=6 on192x384, the maximum interior ratio is4.33; six radial rows and
six polar rows at each axis end have this sign reversal. The opt-in
`HISPID_MONOTONE_PRECONDITIONER=ON` replaces a by max(a,|b|h/2) only in the
finite-difference approximate inverse. This is the minimum diffusion needed
for nonnegative drift/diffusion neighbor weights; it does not assert that
the full potential-bearing matrix is an M-matrix. Spectral residuals, JVPs,
source values, modal basis, stopping norms and physical acceptance are unchanged.
The default remains OFF and the scalar-source experiment is OFF in this trial.

The12x24x16 nonlinear control converges in5 Newton/98 Krylov (baseline6/135).
Residual and JVP witnesses are bitwise identical. Converged sampled fields
differ by1.12e-7 in scaled maximum norm at the existing stopping threshold;
this does not qualify fine-grid equivalence or physical accuracy. Stencil
audit, control and reproduction script are `c4-experiment/monotone-preconditioner`.
Single-GPU shared job59417256 is building this variant for one replay from
the original retained192-grid C4 iterate before any full-solve decision.

Monotone-preconditioner replay59417256 completed in3:10 including build/setup.
The linear solve took59 Krylov iterations and16.365s, reaching true relative
residual7.40468e-4. Step maximum1.42514e12 remains excessive; damping1/1024
is first accepted, with residual ratio0.999441155. This improves the immediate
line-search cutoff by one halving but not the underlying large-update behavior.
No full solve is justified by this replay. All records are in
`monotone-preconditioner/gpu-v1`; image hash
f9ddec38b32fd97519c880ba54f0e0666252f667bf44d26717b137a3c85a0357.
The option stays OFF by default and unqualified for production. No allocation
remains active from this trial.

The next representation question is whether regularity can be imposed as
axis constraints in lower-power modal variables, replacing unresolved
near-axis collocation equations instead of dividing their noisy Fourier
coefficients by higher powers. For the C2 variables with m>=5, requiring
P(t=0,eta)=P(t,eta=+/-1)=0 gives an extra factor t*(1-eta^2) locally and
raises the axis order by two. A tau formulation would require a nonsingular
choice of independent endpoint rows, matching JVP and preconditioner changes,
and manufactured/physical convergence evidence. This is a proposed direction,
not an implemented replacement or an accepted change to physical equations.

A standalone dense flat-scalar m6 experiment now tests this axis-tau proposal
before touching the native solver. Unknowns are nodal C2 P with r=4; polar
endpoint interpolation conditions replace the first/last polar rows for
every radial node, and radial endpoint conditions replace the first radial
row only at remaining polar nodes. This gives2*na+nb-2 independent boundary
rows without duplicate corner conditions. All remaining rows retain the
analytical mapped Laplacian. The physical field is
-2*(1-t)*(sqrt(t)*sqrt(1-eta^2))^6, and its source is checked independently
using Cartesian distance-mode Hessians (maximum scaled discrepancy7.81e-17).

On12x24,20x40,32x64, tau off-grid field errors are5.82e-6,2.93e-9,2.40e-14.
With the same deterministic absolute1e-16 source perturbation, finest-grid
field noise is2.17e-11 versus3.43e-6 for factored C4 (about1.58e5 less).
This supports a native prototype, not binary acceptance: curved coefficients,
vector coupling, nonlinear solve and GPU integration remain untested. The
initial prototype emitted macOS BLAS floating-status warnings despite finite
outputs; explicit einsum contractions and finite checks remove those warnings
and reproduce the conclusion. Both outputs are retained under
`c4-experiment/axis-tau`; `manufactured.json` is authoritative. The standalone
driver is `validation/diagnose_trumpet_axis_tau.py` and uses NumPy only.

`HISPID_AXIS_TAU=ON` now implements the endpoint formulation in the native
residual/JVP with shared CPU/device row algebra and two Kokkos kernels. It
requires C2 modal variables (cap4); m>=5 rows at both polar ends and the
remaining radial-axis edge are replaced by exact Chebyshev endpoint P=0
conditions. Other modes/rows retain their spectral equations. The modal
preconditioner uses matching local two-point endpoint extrapolation rows,
unit scaling and no potential term on those rows. This changes the discrete
boundary formulation, unlike the earlier preconditioner-only experiments.
Defaults remain unchanged. Host/device scratch and optional source-cache
allocation costs are included in memory accounting.

The new identifier is `modal_P_C2tauC4_map_v5_r<lambda>_k<kappa>` and its
residual label explicitly includes axis_tau. Checkpoint validation, map
parsing and initial-guess prolongation recognize this encoding; older native
consumers must reject it. The stored variables retain C2 encoding, so a
verified C2 checkpoint may be used only as an explicitly sourced fresh initial
guess, never inheriting physical acceptance. Unscaled internal diagnostics
now mix PDE and endpoint rows; use the independent physical observer for
physical constraint claims.

The native12x24x16 moderate control converges in5 Newton/69 Krylov. A separate
row-placement check reconstructs the modal endpoint conditions and agrees
to1.11e-16; a centered JVP check has relative L2 error9.74e-12, and a checkpoint
roundtrip preserves coefficients and the exact identifier. Records are in
`axis-tau/native`. Single-GPU shared job59418528 builds the production path,
runs the same structural control on CUDA, then starts the moderate192x384x16
solve only on success. It uses the original source assembly and preconditioner
interior stencil,10 backtracking halvings and unchanged physical bounds.
No native tau physical solution or consumer/horizon acceptance is yet claimed.

The tau CUDA structural control passed in job59418528: endpoint-row error
3.33e-16, directional JVP relative error1.00e-11 and exact checkpoint roundtrip.
Image03abf09a3c2f2100966d21c7462718a47fe918f02725e2017fdc25f1862f942e,
records `axis-tau/gpu-v1`. The moderate fine-grid solve completed setup in
65.148s and is still running; no linear/Newton completion or physical result
is available at this checkpoint.

AthenaK's existing reader already compares the checkpoint identifier exactly
against its linked sampler. Therefore tau integration requires a matching
CPU sampler build, not a weaker reader check. The prepared
`axis-tau/build_consumer_perlmutter.sh` builds that sampler and the existing
AthenaK pgen in separate output directories. It has not been launched.

The first native tau fine-grid run59418528 completed8:32, but its first linear
solve exhausted2400 iterations at relative residual0.008948 versus0.001; no
Newton update was accepted. Solve235.606s, workflow406.938s. The retained
zero-correction checkpoint is7c01f0713fd4e3d32d0d881ae7fb6a3c09990c0b417053a01d99dd76ee7f633b.
Near H/M RMS0.1874/0.04975 and bulk0.000941/0.001329 fail physical bounds.
Full logs and result are `axis-tau/gpu-v1`; no horizon run is warranted yet.

The radial tau rows of the preconditioner now use exact endpoint weights:
its existing dense radial blocks can store these with no additional matrix
size. Polar endpoint rows still use local two-point approximations to retain
block-tridiagonal storage. Residual and JVP witnesses remain bitwise equal
to the original tau path. The small nonlinear control takes5 Newton/61 Krylov
versus5/69, with converged sampled fields agreeing to1.01e-13. An isolated
CUDA build/replay from the zero-correction checkpoint tests the fine-grid
linear solve before any further full solve. Records and command are
`axis-tau/radial-exact`; no physical acceptance is claimed.


Exact-radial replay59419101 completed6:40 but exhausted2400 Krylov steps at
true relative residual0.01616844 (target0.001), with no Newton update;
solve235.035s. Image6019e93f6af829e437854a46446daab6ab789a855f645d8a9c8abca97d18582b.
The JSON, trace and image hash are retained in `axis-tau/radial-exact`.

`HISPID_EXACT_POLAR_TAU` (off by default, requires axis tau) now adds exact
polar endpoint rows through a rank-2na Woodbury update. `HiSpID_modal_block.hpp`
shares the existing base factors and batches their matrix-RHS solves; CPU
and Kokkos apply the same response/Schur correction. Only the preconditioner
changes; defaults, physical residual/JVP and checkpoint representation remain.
Memory guards account for response banks, Schur setup and device copies.
The independent dense-matrix test uses Chebyshev-series endpoint evaluation,
pivoted factors, nonzero boundary data and interleaved component isolation:
solution error<=6.67e-15, scaled matrix residual<=1.85e-15. The small coupled
CPU solve takes5 Newton/53 Krylov, versus5/61 for exact radial alone;
residual/JVP witnesses are bitwise unchanged and sampled field/gradient
relative difference2.53e-13. Records: `axis-tau/polar-exact`.
Perlmutter shared-GPU job59420201 runs the matrix and small coupled GPU controls,
then only on success the original fine-grid zero-state replay. No full binary
or physical-convergence claim is made from these controls.

The polar-correction CUDA control passed in job59420201:5 Newton/54 Krylov,
field/gradient agreement2.82e-13 with the CPU result. Matrix control also
passed. Image9489bca44c0bbd6ea20ae4996b4e2cdf225e84e3436ec80b253c4ffaf2b0c3e6;
records `axis-tau/polar-exact/gpu-v1`. The job has entered the fine-grid replay;
its terminal result is not yet available. Implementation commitede4316 pushed.

Fine-grid replay59420201 completed3:13. The exact polar correction reaches
the0.001 linear target in6 Krylov steps, true relative5.73775e-4. Full Newton
step accepted, residual ratio0.02631594, positive minpsi1.00000668; solve13.008s.
The one-Newton-step diagnostic intentionally reports its iteration limit;
this does not indicate a failed linear solve or establish physical acceptance.
Raw trace/JSON retained in `axis-tau/polar-exact/gpu-v1`.
This justifies full moderate192x384x16 zero-start solve59420346, same pinned
image, forcing0.001 and unchanged physical bounds. Command:
`axis-tau/polar-exact/full_perlmutter.sh`; remote output
`axis-tau-polar-full-v1/moderate192`. Allocation is one shared GPU,20minutes.

Full moderate192 tau run59420346 completed4:15. Nonlinear convergence in
5 Newton/53 Krylov, all full steps; weighted Linf max5.49e-14. Solve65.395s,
setup65.537s, workflow246.599s, peakRSS6317476KiB. Independent near H/M RMS
1.128066e-6/3.789400e-7 and bulk5.738737e-8/4.547640e-7; near Hamiltonian
alone exceeds1e-6, so this grid remains diagnostic. Positive sampled psi
and metric, exterior stencils unmodified. Checkpoint
61b7eb4d2fe0aed07c97957193aa7fca25722a9b82fce3c4e374be2d7c55c98f;
records `axis-tau/polar-exact/full-v1` (large coefficients retained remotely).

Single-GPU shared sequence59420603 now runs224x448,240x480,256x512, nphi16,
same image and physics, using the preceding converged state only as an
initial guess. Restart24 avoids exceeding the unchanged64GiB aggregate
allocation guard with the new border matrices. Necessary resolution checks
use one production pipeline, not a backend matrix; each level must pass
physical bounds before proceeding. Final three-grid assessment checks raw
observers, charge stability and decreasing RMS. Script:
`axis-tau/polar-exact/sequence_perlmutter.sh`; remote `axis-tau-polar-sequence-v1`.

CPU job59420409 completed3:26: matching sampler and AthenaK consumer built.
Build images/hashes are retained in `axis-tau/polar-exact/consumer-build`;
script `build_polar_consumer.sh`. New tau-checkpoint import/horizons have not
yet been run. Condensed obsolete failed-trial prose in the LaTeX report into
a comparison table; all historical raw results and this ledger remain intact.

Prepared and staged `axis-tau/polar-exact/horizons_perlmutter.sh` for the new
finest tau checkpoint. It verifies the pinned GPU producer, CPU sampler and
AthenaK images, checks completed per-grid physical bounds, runs independent
sampler migration, then strict l8/12/16 plus quadrature horizon checks with
16 geometry threads and factorized harmonic storage. It has not been
launched: sequence59420603 is still running its224 level. This avoids a
redundant consumer/horizon campaign on each intermediate grid. Remote script:
`/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005/horizons-polar-tau.sh`.

Sequence59420603 failed10:36 at the224 gate, and did not launch240/256.
It accepted two undamped Newton steps, then exhausted2400 Krylov iterations
on the third linear solve (true relative0.0244754 vs0.001). Total3 attempted
Newton/2428 Krylov, nonlinear max4.33725e-12>1e-12. Nevertheless, all independent
physical bounds pass: near H/M RMS2.08042e-7/1.68057e-7, bulk4.87446e-8/5.03784e-7.
This is retained as a failed solve, not accepted by its physical checks alone.
Solve392.298s, workflow625.940s, RSS9930348KiB. Checkpoint
7687a64d90a29e1f2272255e99c91a6bde01f7c44aa3ab51cd858861bbc17721.
Records `axis-tau/polar-exact/sequence-v1` preserve JSON, physical observers,
trace and job id; large checkpoint/coefficients remain remote.

Job59421025 resumes that retained224 iterate with the existing standard
linear forcing0.1. This targets the observed inner-solve stagnation without
changing final nonlinear tolerance1e-12, native equations/image, free data,
resolution or any physical/horizon acceptance criterion. It continues to240
and256 only after completed nonlinear and physical gates pass. New outputs
`axis-tau-polar-sequence-v2`; script `sequence_standard_forcing.sh`. The prepared
horizon workflow now points to this sequence's finest checkpoint and has not
been launched. No duplicate of the still-running old job was created: its
terminal FAILED state and released allocation were verified before restart.

The standard-forcing continuation224 in job59421025 completed its level and
passed nonlinear and physical gates. One Newton/119 Krylov, solve36.829s;
weighted max3.34966e-13<1e-12. Near H/M RMS2.07100e-7/1.65066e-7,
bulk5.02305e-8/4.33618e-7. Workflow266.592s; no physical threshold changed.
Checkpointd78a1b3f8cc8d517f9ef33b88258981c3669ab438647a0c3d0ea2353b98ca735.
Records `axis-tau/polar-exact/sequence-v2/moderate224`.
The job has advanced to240. This is one qualifying grid, not completed
three-grid convergence or horizon acceptance.

While240 remains live in59421025, checked the existing allocation bound with
actual Cached size1984 bytes (all-double structure, stable scalar sourceOFF).
At240x480x16, restart64 estimates60.71624GiB; at256x512x16, restart32 estimates
63.55018GiB, while restart64 would exceed64GiB. This corrects the deliberately
overestimated4096-byte cache used in the initial preflight; no guard or native
code is changed. `restart-memory-preflight.json` records source hashes and
estimates, not measured peak memory. A larger restart is a feasible targeted
retry if the live240 linear solve fails; it has not been applied to that run.

Sequence-v2 job59421025 terminated after failing240;256 was not launched.
The240 level accepted four full Newton steps, then exhausted2400 Krylov steps
at relative0.8729704 vs0.1. Total5 attempted Newton/4311 Krylov, nonlinear
max1.61556e-9>1e-12. Near H/M RMS7.51953e-7/2.19342e-7; bulk1.12602e-7/1.72401e-6
also fails momentum bound. Solve824.584s, workflow1084.748s, RSS12152192KiB.
Checkpoint8c920fe6ef5ed7392db21e6cbea8df5a064b7a47aad2a137f80c70085810a275.
Records retained under `axis-tau/polar-exact/sequence-v2/moderate240`.

Job59421758 resumes that state with restart64 at240, keeping forcing0.1 and
all final criteria unchanged. Existing allocation estimate60.71624GiB fits
64GiB; if240 passes,256 uses restart32 (63.55018GiB). No guard modification.
The passed224 result is reused by the final three-grid assessment. Outputs:
`axis-tau-polar-sequence-v3`; command `sequence_larger_restart.sh`. Old job
terminal failure/released allocation verified before this single-GPU launch.
Prepared, unlaunched horizon script now targets this sequence's256 dataset.

Job59421758 terminated FAILED11:40 and released its allocation. Restart64
reduced the failed240 inner relative residual to0.1210990354, still above0.1
after2400 iterations; no Newton update was accepted. Nonlinear max1.61558e-9,
bulk momentum RMS1.72401e-6 remain failures. Solve430.714s, total692.763s,
RSS12098548KiB. No256 run launched. Checkpoint
 d852ca6013190b37fb395c342eb0559767d6b1adba0ab21d83cc387fe45066fc.
Raw result/physical observers/log preserved in polar-exact/sequence-v3/moderate240;
large state artifacts remain remote. The next step is targeted diagnosis of
linear operator/preconditioner error; no further blind restart/forcing sweep.
Source review confirms GMRES uses two-pass orthogonalization and explicitly
recomputes the true residual at restarts and termination. The current Newton
trace starts only after a successful inner solve, explaining the empty trace
while this failed solve was running; it does not imply lack of computation.

Added opt-in HISPID_PROBE_MODAL_FACTORS=1 diagnosis. Each assembled modal
preconditioner captures its original nonzero matrix entries before factorization,
forms a deterministic all-frequency manufactured RHS with long-double accumulation,
and reports inverse forward error, componentwise backward error and relative L2
residual. Exact polar endpoint rows are included in the target matrix. Only O(N)
sparse diagnostic storage is added; the production default is unchanged. Group
2m is scalar and2m+1 is the shared vector factor. The final Krylov recurrence and
explicit true residual are also printed, permitting a recurrence-gap diagnosis.
Small reference residual/JVP/fields/linear-history remain bitwise identical;
existing independent dense polar-border controls pass with the probe enabled.

Job59422127 runs one bounded64-Krylov replay of the failed240 checkpoint, using
a separately built diagnostic image/source directory and one shared GPU. It does
not run physical observers or launch a new resolution. Records/scripts are in
polar-exact/factor-probe; remote outputs axis-tau-factor-probe-v1. No result yet.

Factor-probe job59422127 completed4:15, allocation released. All18 host modal
factors recover the manufactured witness with forward error<=1.11745e-9
(worst scalar m6); componentwise backward error<=3.17190e-10 and relative
L2 residual<=7.27731e-16. These are witness-specific accuracy measurements,
not condition-number bounds. After64 GPU Krylov steps, true relative residual
is0.92150255; recurrence2.9788492841399357e-8 and explicittrue
2.9788492841393478e-8 differ relatively1.974e-13. This is real slow convergence,
not a recurrence-only false residual decrease. GPU application of the inverse
has not been isolated from the approximation quality and remains the next
specific diagnostic distinction. Native image9b3eb7bb24020bc497095223e30ff748d0b1b64c832f88ca3169d64fd552b3b6.
Raw records and parsed per-group analysis in factor-probe/gpu-v1. No physical
acceptance is claimed from this deliberately iteration-capped diagnostic.

Added opt-in HISPID_PROBE_LINEAR_ACTION=1 to distinguish device inverse error
from preconditioner approximation. It compares the actual Newton RHS through
host triangular/Woodbury solves and GPU inverse/Woodbury kernels, reporting
absolute differences per Fourier mode/component. After the bounded Krylov run,
it resolves the initial and remaining true residual norms by Fourier mode and
component, including the tau-row share. The additional JVP is counted as
work; neither the operator nor the solver tolerance is changed. Small reference
residual/JVP/field witnesses and iteration histories remain bitwise equal.
Job59422392 is the bounded64-step240 replay with these diagnostics on one shared
GPU, using separate source/build directories axis-tau-action-probe and output
axis-tau-action-probe-v1. No new physical resolution or acceptance claim.

Action-probe job59422392 completed4:08 and released its allocation. The
actual-RHS GPU/host preconditioner comparison differs7.48324 percent globally
in L2 (worst mode/component7.94297 percent), dominated by mode4 transverse
momentum. High tau modes alone agree closely; looking only at those would
misdiagnose the result. Reference mode4 transverse correction norms are227.0
and227.1, with device differences15.88 and18.04. Remaining true-residual energy
is99.7173 percent momentum and70.9224 percent modes0/1; tau rows contribute
only7.002e-16 of energy. The next distinction is Fourier projection/row scaling
versus CUDA inverse application on identical projected input. Do not replace
the mathematical preconditioner before isolating this concrete discrepancy.
Raw records and parsed analysis in action-probe/gpu-v1. No full solve launched.

Factored the existing host Fourier/row-scaling operation into project_rhs,
without changing its arithmetic. Modal::apply now optionally accepts already
projected/scaled input (defaultfalse, all production callers unchanged). The
opt-in action probe compares both ordinary device inversion and inversion of
the exact host-projected RHS against the same host solution. This isolates
projection from inverse-kernel differences without duplicating either inverse.
Small reference residual/JVP/fields/iteration histories remain bitwise equal.
Job59422608 uses one shared GPU and only1 Krylov step; its useful observations
are the two pre-solve actual-RHS comparisons. Separate source/build directories
axis-tau-matched-probe preserve the earlier images. Outputs are
axis-tau-matched-probe-v1; no numerical result yet.

Matched-projection job59422608 completed4:04, allocation released. Ordinary
GPU/host inverse relative difference0.0748324 falls to2.98125e-16 when both
consume the identical host-projected/scaled RHS. This isolates Fourier
projection, not inverse application. Raw records in matched-probe/gpu-v1.

Added opt-in HISPID_COMPENSATED_MODAL_PROJECTION (defaultOFF). A shared
host/device helper compensates mean subtraction, product error via FMA, and
summation; HiSpID selects it for its modal preconditioner. The shared device
Modal default remains uncompensated, preserving BY behavior. Equations and
checkpoint basis are unchanged. Exact dyadic cancellation tests verify lost
product and subtraction bits and constant annihilation; host controls pass.
Small nonlinear control5Newton53Krylov, residual/JVP bitwise equal, sampled
fields scaled difference3.11673e-13. Job59422875 builds the GPU implementation,
runs the same exact cancellation tests on-device, then resumes the saved240
state with restart64/forcing0.1 and all final criteria unchanged. Actual-RHS
host/device diagnostics remain enabled during this corrective replay. Remote
source/build axis-tau-compensated, outputs axis-tau-compensated-v1; local
controls/scripts in polar-exact/compensated-projection. No GPU result yet.

Job59422875 built successfully and passed the three exact cancellation
witnesses on both CPU and CUDA (answers -2^-54,2^-53,0). Pinned image/source
hashes and cancellation output are retained in compensated-projection/gpu-v1.
The resumed240 physical workflow remains live; no solve result yet.

The first corrected replay59422875 failed3:46 before solving: its command
omitted --initial-source-library, so the provenance guard correctly refused
the old-image checkpoint. Retained failed-launch.log and incomplete-result.json;
no numerical failure or accepted step occurred. Job59422941 reuses the built
and qualified image7398b4b254371c39409d511ebe0808fef5e8c45b9c1e6b1af2ccb52acecf31ea,
explicitly pins old source image9489bca4..., and passes its path via the existing
source-library interface. Same maps/grid use the identity prolongation path;
all physical free data and final criteria are unchanged. Outputs
axis-tau-compensated-v2; script resume_perlmutter.sh. Previous allocation release
verified before the single-GPU retry. No rebuild or duplicate test campaign.

The corrected first actual-RHS comparison in live job59422941 gives global
GPU/host inverse relative L2 difference2.816881431e-16, down from
0.0748324. Worst mode/component relative difference3.918212561e-11.
Initial36-row snapshot preserved in compensated-projection/gpu-v2. This
qualifies the targeted projection correction on the failing240 state, not the
nonlinear solve or physical acceptance, both still pending.

Corrected replay59422941 finished13:12 and released its allocation. The job
completed its diagnostics, but the solve did not converge: first inner solve
319 iterations/relative0.0917176 accepted a full Newton step; second inner
solve2400/relative0.305647 failed target0.1. Total2attemptedNewton/2719Krylov,
solve524.666s, workflow787.662s, RSS12275884KiB. Nonlinear max2.16684e-11
still exceeds1e-12. Near H/M RMS5.23509e-7/2.57229e-7; bulk5.40303e-8/2.53470e-6
still fails momentum bound. Saved checkpoint
2a354f2337771e9c4e40d5750fae0a880d64513379c789be22a091fe158c4042.
Final JSON/raw physical observers/log retained in compensated-projection/gpu-v2.

Job59423304 resumes this improved state with inexact-Newton forcing0.5,
restart64 and the same image/grid/physical free data. This is justified by the
measured true linear residual0.306, not a recurrence-only estimate: if
Js=-r+e and ||e||<=eta||r|| then r.Js<=-(1-eta)||r||^2. Eta0.5 admits a descent
direction, and the existing line search verifies actual nonlinear decrease.
The final nonlinear1e-12 and physical/horizon thresholds are unchanged; no
failed solution is accepted. Routine inverse probes are disabled after the
actual-state correction was qualified. Script inexact_newton_perlmutter.sh;
remote output axis-tau-compensated-v3. A final explicit per-grid gate now
makes job failure reflect nonconvergence or failed physical bounds. No rebuild
or additional resolution launched.

Continuation59423304 finished its diagnostic workflow in813.315s but failed
nonlinear and physical acceptance. Forcing0.5 admitted one full step after569
Krylov iterations (true relative0.499639); the second inner solve exhausted2400
at0.546869. Total2969Krylov,552.286s solve; nonlinear maximum2.13643e-11
still exceeds1e-12. Near H/M RMS5.21125e-7/2.57128e-7 pass, but bulk momentum
RMS2.53426e-6 fails. No threshold is relaxed and no finer grid is launched.
Raw physical observers, full log and final JSON are retained under
polar-exact/compensated-projection/gpu-v3. Saved diagnostic checkpoint SHA256
 e997e7af3cc437b39d2d533a208163666a53d5f7bbcb07e94127db155de6666a.
The relaxed forcing does not cure linear stagnation; further work must address
that difficulty rather than treating another tolerance relaxation as validation.

Job59423677 runs a bounded one-Newton BiCGStab probe from the final v3
checkpoint on the identical corrected image. The existing shared BiCGStab
backend is selected explicitly by a new diagnostic-driver argument; default
GMRES is unchanged. Limit1200 BiCGStab iterations allows approximately2400
JVPs, comparable to the failed GMRES step, with target0.1. This isolates a
potential restart-related difficulty without another full physical campaign;
it cannot by itself prove that restart is the cause. No equation/preconditioner
change or new library build. One shared GPU,30-minute allocation; output
axis-tau-bicgstab-probe-v1, script bicgstab_probe_perlmutter.sh. Driver help and
shell syntax checked. Physical acceptance remains false for diagnostic output.

BiCGStab probe59423677 completed but did not produce a usable Newton step:
1200 iterations/2402JVPs/2400preconditioner applications,401.789s solve,
true relative residual1.73413e103. No update was accepted. This rejects the
existing BiCGStab backend for this particular frozen240 system; it does not
establish the cause of GMRES stagnation. Raw probe JSON/log/job ID retained
in compensated-projection/bicgstab-probe. Next investigate retaining useful
GMRES directions across restarts while preserving true-residual checks and
the memory bound; do not repeat the backend matrix or relax physical gates.

An opt-in shared PK_LGMRES controller now appends up to three normalized
previous-cycle corrections to each ordinary restart subspace, recomputing
operator images and retaining true-residual verification. Source reference:
https://docs.scipy.org/doc/scipy-0.18.0/reference/generated/generated/scipy.sparse.linalg.lgmres.html
This is right-preconditioned augmentation of the existing shared Arnoldi
controller, not a SciPy dependency or cross-Newton recycling. Existing method
values, option/result struct layouts and default arithmetic are unchanged.
Three saved directions plus a correction scratch and six extra Arnoldi/search
vectors require at most ten additional vectors. HiSpID selection/memory guards
and CUDA qualification are not yet wired, so this is not a production change.
CPU independent controls cover all three methods, failures, zero RHS and
initial guesses:227 checks pass, also under AddressSanitizer/UBSan. A separate
nine-dimensional diagonal spectrum(.01,.02,.04,1,2,4,8,16,32), exact solution
all ones, gives ordinary GMRES(5) residual5.53915e-7 after300 iterations;
augmented GMRES with five new directions plus up to three saved directions
converges in90 iterations to true4.26767e-13. This is an algorithm control,
not an equal-memory performance claim or evidence about the binary system.
A direct C-as-C++ compile attempt failed on existing C void-pointer conversions;
it is not a Kokkos compilation test. CUDA qualification remains pending.

HiSpID now explicitly selects lgmres through native/Python options, preserving
old enum values/ABI and default choices. Its allocation guard adds ten vectors
plus1MiB scalar headroom to the existing aggregate bound before solving.
A20x28x8/restart64 context fits40MiB but the augmented solve is correctly
rejected. The production240/restart64 estimate is61.266533GiB below64GiB.
Small CPU nonlinear restart5 augmentation:5Newton68Krylov, sampled field
scaled difference9.46951e-14 from the retained compensated baseline; residual
and JVP witnesses bitwise unchanged. Records in polar-exact/augmented-gmres.
Job59424143 builds the selected CUDA implementation in source/build-axis-tau-lgmres,
runs only shared Krylov controls and the small augmented HiSpID solve, then
probes one Newton step of the saved240 state with restart64/target0.1/max2400.
No physical acceptance is inherited. Script run_perlmutter.sh; output
axis-tau-lgmres-v1. CUDA results pending. A local follow-up guard explicitly
rejects LGMRES in the BY-specific Kokkos wrapper (its capacity formula supports
only the existing two methods); this guard is not in the already dispatched
source snapshot and is irrelevant to the HiSpID-only probe. It will be included
in the next source build. BY's Python adapter continues to reject lgmres.

CUDA qualification in live job59424143 passes221 shared Krylov checks. The
separated-spectrum witness converges in90 iterations, independent residual
4.32480e-13. Small augmented HiSpID solve:5Newton68Krylov,0.07273s solve,
CPU/GPU sampled-field scaled difference2.57439e-13. Pinned CUDA image
2cfec4d71e54542add9805ee544dc73fbe4f54c3af8948e8ba1aabfabc61cdf1.
Control records retained in augmented-gmres/gpu-v1. The one-step240 probe
has started; no fine-system convergence or physical acceptance claim yet.

The full240 augmented probe59424143 stops after2400 iterations at true
relative0.47081893 (target0.1),421.3998s solve,2437JVPs/2301preconditioner
applications. No Newton update was accepted. This is modest improvement over
the same-state GMRES residual0.546869 but does not resolve stagnation. Final
probe JSON/log retained alongside qualification records. Further Krylov-method
sweeps are not justified; inspect the actual operator/preconditioner and spatial
error structure before another expensive solve. Physical acceptance remains
unchanged and false.

The centered-stencil audit at224/240 finds no negative interior neighbor
weights for modes0/1, and the same negative-row pattern for mode4 at both
resolutions. This does not explain the resolution-specific failure; do not
launch a monotone-preconditioner trial solely from that sign pattern.
Reanalysis of retained corrected-GMRES diagnostics (v2 final failed step)
shows98.53% momentum residual energy, but only26.57% in modes0/1, unlike the
older pre-correction70.92% result. Remaining energy spans modes0--7;
tau-row fraction1.94e-13. Records/scripts in polar-exact/residual-localization.
Job59424678 evaluates the saved224 and240 residual once each with the same
qualified CUDA image2cfec4d7..., no Newton or Krylov solve. It records spatial
energy concentration, component fractions and endpoint layers. Checkpoint
source hashes and basis identities are checked; physical acceptance is not
inherited. Script run_perlmutter.sh, remote axis-tau-residual-localization-v1.
The observations will distinguish endpoint concentration from bulk operator
error before changing the preconditioner. No new build or resolution sweep.

Localization59424678 completed3:25 and released its allocation. The same CUDA
image evaluates224 with L2=3.54850e-11/Linf=3.34966e-13 and240 with
L2=1.48138e-9/Linf=2.13643e-11. Momentum carries98.90%/98.88% of energy.
For240, the outermost eight radial/polar index layers contain only4.99e-6/
3.29e-8 of residual energy. Largest point is(i,j,k)=(51,55,4); energy is
spread across interior indices and azimuths. This does not support endpoint
rows as the sole remaining problem. Raw224/240 JSON and job ID retained in
residual-localization/gpu-v1. Before changing the preconditioner, compare
actual-state residual evaluation precision between CPU and CUDA: a fine-grid
roundoff floor has not yet been excluded by the small-grid controls.

Job59424759 compares the actual failed240 residual between CUDA and native
CPU reference execution in sequential contexts of the same qualified image.
No solve, tolerance change or preconditioner modification. This tests whether
backend evaluation differences are significant relative to the retained
1.48e-9 L2 residual, rather than extrapolating small-grid agreement. Neither
path is treated as exact truth. The driver retains both raw residual arrays
with a hash for later analysis, along with component differences and unknown
magnitudes. Same source checkpoint hash and parameterization guards; output
axis-tau-residual-reference-v1, script reference_perlmutter.sh. Previous
localization allocation release was verified. Result pending.

At10:11, job59424759 remains RUNNING. Code inspection confirms that native
reference context creation builds geometry serially, unlike the production
host Kokkos path. An attempted extension from15 to35minutes was denied by
Slurm; the live job was left intact. The diagnostic driver now saves a hashed
CUDA-only partial result before starting the reference context and emits
flushed stage progress. This change is for subsequent invocations; it does
not alter the running process or any numerical implementation. Partial files
explicitly indicate reference_pending and cannot establish backend agreement.
Python syntax and diff checks pass; no numerical rerun was added for logging.

Job59424759 terminated TIMEOUT at15:10; session21810 exited143 and no final
comparison artifact was produced. Its allocation was released before the
replacement job59425181 started (session31839). The replacement uses the same
image/checkpoint and diagnostic, with45minutes requested and stage-saving
driver ae0e037. Output axis-tau-residual-reference-v2; the script accepts
PUNCTURE_RESIDUAL_OUTPUT to preserve the first attempt. No solve or numerical
parameter was changed. This repeats the unfinished diagnostic because of
insufficient wall time, not an additional backend/performance sweep.

Reference comparison59425181 COMPLETED27:29; session31839 exited0 and the
allocation was released. CUDA setup/residual completed at99.20s; serial
reference context became ready at1616.01s. Same checkpoint/image residuals:
CUDA L2=1.481382400533222e-9, CPU L2=1.4813825027858145e-9;
difference L2=2.1104656688840964e-13, Linf=2.5677291872716226e-15,
relative to CUDA=1.4246596072185256e-4 (0.01425%). This does not support a
dominant backend-specific evaluation error; agreement does not exclude
common discretization/precision errors. No additional backend comparison is
justified by these results. Next inspect the operator/preconditioner mismatch
rather than another Krylov or tolerance sweep. Physical acceptance unchanged.
Raw JSON/log/job ID are in residual-localization/reference-v2. Both arrays
remain in remote axis-tau-residual-reference-v2/residual240.npz with SHA256
1f52fef6a78aa4d8c9867c74ab99239474365c5089072aa8d2aad8c41fe8d414.

Reanalysis of those saved arrays (no new solve) finds a grid-scale residual:
orthonormal radial/polar DCT-II median degrees238/477 on240/480 nodes.
The last eight radial/polar modes carry93.59864%/78.55749% of squared residual;
upper halves carry98.66988%/87.67540%. CPU reference gives the same pattern.
Script residual-localization/spectrum.py validates the array hash and Parseval
identity (maximum discrepancy2.0e-15); a known product of degree3 radial and
degree7 polar cosine modes is recovered to1e-13. Result reference-v2/spectrum240.json.
These are weighted residual spectra, not spectra of the solution or independent
physical constraints. They identify grid-scale stagnation but do not prove its
cause. A targeted next candidate is exact spectral radial differentiation in
the already-dense radial preconditioner blocks, retaining polar FD structure,
tau conditions, memory guards, and the unchanged residual/JVP. This tests
high-frequency approximation quality without filtering or another method sweep.

Opt-in HISPID_SPECTRAL_RADIAL_PRECONDITIONER now replaces the radial FD
derivative contribution with exact polynomial differentiation transformed
through t=lambda*s/(1-(1-lambda)*s), s=(1+z)/2. Off-diagonal differences
annihilate constants. Polar FD and tau rows are unchanged; scalar/vector
factors use the existing dense-radial block storage on CPU and CUDA.
CPU uses compact assembly for this option to avoid dense rows in per-point
sparse containers. Added host derivative/map tables cost16*na^2+24*na bytes
and are included in the aggregate guard. DefaultOFF; no BY change.
The constant/degree3/degree(n-1) manufactured operator control over12/32 nodes,
three radial stretches and regularity exponents0/1/4 has worst scaled
error7.00031e-13. Small native binary converges5Newton48Krylov (baseline53);
sampled-field difference1.91263e-13 and residual/JVP witnesses bitwise unchanged.
Caching radial map factors retains identical sampled fields and iterations.
Records and production qualification/probe script are in polar-exact/spectral-radial.
Fine-grid improvement is not yet established. The planned single-GPU run
checks the changed small solve against this CPU result before a frozen240
GMRES probe with unchanged restart64, target0.1 and budget2400.

Job59426128 (session43880) is building the spectral-radial image on one shared
80GiB GPU allocation with25minutes requested. Remote source/build suffix
axis-tau-spectral-radial, output axis-tau-spectral-radial-v1. The small CUDA
solve must converge and agree with retained CPU fields within1e-10 before
the frozen240 probe executes. No other compute allocation remains active.
