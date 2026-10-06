# Trumpet implementation and acceptance ledger

Goal started 2026-10-05. Native base4ce5974; AthenaK consumer base9495fdb9.
Dedicated branches: codex/hispid-trumpet and codex/hispid-trumpet-pgen.
The prior worktrees and raw results are unchanged. No previous benchmark
campaign is resumed. Current status: seed-family API, context/cache/sampler
dispatch and versioned checkpoint integration implemented; CUDA controls and
full consumer build pass. Isolated spin.99 horizon convergence passes; two moderate binary grids
converge internally but fail physical accuracy. No binary acceptance yet.

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
| Moderate binary | Three increasing resolutions; independent fixed exterior near/bulk H/M RMS<1e-6 and max<1e-4 with decreasing errors; positivepsi; charge changes<1e-3 | one passing grid; convergence pending |
| Consumer/horizons | Bound checkpoint and field/gradient roundtrip<1e-12; initial-time expansion RMS<1e-7 at3 angular orders, mass/spin changes<1e-4; modified-region enclosure | isolated spin.99 passes; binary pending |
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
