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
