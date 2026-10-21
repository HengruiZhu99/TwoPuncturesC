# Trumpet implementation and acceptance ledger

Goal started 2026-10-05. Native base4ce5974; AthenaK consumer base9495fdb9.
Dedicated branches: codex/hispid-trumpet and codex/hispid-trumpet-pgen.
The prior worktrees and raw results are unchanged. No previous benchmark
campaign is resumed. Current status: seed-family API, context/cache/sampler
dispatch and versioned checkpoint integration implemented; CUDA controls and
full consumer build pass. Isolated spin.99 and the finest moderate binary
horizon checks pass, including angular convergence and modified-region
enclosure. The retained three-grid QI sequence meets the physical error bounds,
but fails decreasing bulk errors. No full binary acceptance yet.

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
| Moderate binary | Three increasing resolutions; independent fixed exterior near/bulk H/M RMS<1e-6 and max<1e-4 with decreasing errors; positivepsi; charge changes<1e-3 | retained QI grids meet bounds; strict bulk convergence fails |
| Consumer/horizons | Bound checkpoint and field/gradient roundtrip<1e-12; initial-time expansion RMS<1e-7 at3 angular orders, mass/spin changes<1e-4; modified-region enclosure | isolated spin.99 and finest moderate binary pass |
| QI compatibility | Original batch-seed API; background/cache/sampler numerical equivalence is separately assessed after stable arithmetic | seed exports bitwise preserved; historical sampler fixture not transferred |
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

The explicit family API now selects QI=0, R0=M trumpet=1 or Liu=2 for both active
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

GPU qualification in59426128 passes: mapped-polynomial error3.65731e-14;
small nonlinear solve5Newton48Krylov, CPU/GPU field difference2.39090e-13.
Pinned GPU image2dfa13fc56247439e8037ac0b4a857e4821aeabb9c39bf853fa6b9514ecb08b8.
Records saved in spectral-radial/gpu-v1. The frozen240 probe is running;
no fine-grid convergence or physical acceptance is inferred from these controls.

## Human-requested wrap-up — 6 October 2026

Job 59426128 finished the frozen240 spectral-radial probe: failed true relative
residual 0.7866735212674494 against 0.1 after 2400 iterations, 2439 JVPs,
448.510415664 seconds. The qualified small controls do not establish fine-grid
effectiveness. No accepted Newton update or binary acceptance follows.
The probe JSON is retained in spectral-radial/gpu-v1/probe.json.
No additional experiment is launched during this wrap-up. The report is updated
to APS/PRD REVTeX preprint format and exported as the current PDF.
The overall physical-validation goal remains incomplete; earlier results and
the discontinued performance campaign are preserved.

## Moderate convergence repair investigation — 18 October 2026

Resumed at the user's request; the scope is one fixed moderate formulation,
three resolutions, then the finest accepted solution in AthenaK. No extremal
parameter study or performance matrix is reopened.

Candidate `moderate-wide` retains the original 0.15 isolated-radius flat plateau
and widens the transition outer radius from 0.5 to 0.85. This is an explicit
change to the interior extension, not a claimed preservation of the entire
previous solution. The exterior vacuum CTT equations and all physical bounds
are unchanged. Support radii are 0.4684236866 and 0.3143675559; previous moderate
horizons only screen these as plausible. Fresh measured enclosure is required.
A two-context frozen-field control finds exterior residual/JVP bitwise identity
on 15408 collocation points, with changes on 2224 interior points.

The already qualified compensated-projection, exact-polar-tau GPU image
7398b4b254371c39409d511ebe0808fef5e8c45b9c1e6b1af2ccb52acecf31ea is
pinned; the unsuccessful spectral-radial option is not used. A cold160 pilot
runs in job59952590. Job59952530 stopped before creating the solver because
the first workflow bundle omitted checkpoints.py; the corrected bundle includes
all validation Python dependencies and uses a fresh output directory.

The declared candidate sequence is160x320x16,192x384x16,224x448x16 with
unchanged forcing0.1, restart64, nonlinear1e-12 and physical RMS/max1e-6/1e-4.
The existing independent observer and strict three-grid assessor are reused.
Multigrid is deferred: retained stalled residual energy is concentrated at the
highest meridional frequencies, whereas the existing FD approximate inverse
already solves its global operator. No evidence yet justifies an additional
coarse correction. Evidence and gated solve/consumer scripts are in
validation/trumpet/moderate-wide/. This candidate is not yet accepted.

The160 pilot completed:10Newton/47Krylov,58.7603 solve seconds, maximum
weighted residual2.73455e-13. The legacy observer's .6m shell lies inside
the widened modified region; its recorded vacuum values are not an exterior
test. It is retained unchanged. The declared transition-aware profile moves
only such fixed shells to max(legacy radius,1.1*inner_max). It keeps48near/18bulk
points, the original1.5/3 isolated-radius shells and bulk coordinates, all
physical tolerances and the unmodified-stencil requirement. This profile is
frozen across the new three-grid sequence; its near norms must not be presented
as the same spatial norm as the original narrow-transition profile.

A sampling-only reobservation of the identical160 checkpoint took28.75seconds
and passes all exterior bounds: nearH/M RMS7.64989e-7/3.40705e-7, bulk
3.01503e-8/1.49311e-7. The checkpoint hash is
18d79e54c8bc5276006efec25a516ef8c23d3861ec1115433a2f09665ce81630.
Job59953429 exited before sampling because the reobserver lacked the examples
import path; fixed and verified locally before relaunch. Job59953525 runs the
passed reobserver and192/224 sequence on one shared80GiB GPU. The160 solve is
not repeated; coarse metadata refers directly to the original checkpoint.

The192 continuation passes all per-grid bounds:5Newton/34Krylov,64.9588solve
seconds,253.594total,6.11GiB peak host RSS. NearH/M RMS are
1.76064e-7/1.65403e-7 and bulk2.36529e-8/5.28148e-8; all four decrease
from160. The224 solve converged in6Newton/253Krylov (156.094solve seconds,
maximum weighted residual4.49985e-13); its physical observer is pending.
Before receiving that observer, the final declared study was extended to
192/224/240 so it tests the original stalled240 resolution directly. The160
level remains a pilot. All four use identical free data, native images, forcing,
restart, nonlinear tolerance and transition-aware observer rule. The new240
script and final consumer target are recorded; no new allocation for240 has
yet been launched. The consumer uses the previously qualified shape-reuse
option only as an initial guess; every surface/quad/enclosure check is fresh.

The160/192/224 assessment FAILS despite all per-grid bounds passing. At224,
nearH/M RMS4.76322e-8/1.19652e-7, bulk3.20071e-8/9.22180e-8; both bulk
values rise from192. Charge changes remain below8.1e-7. No240 run or horizon
run is launched. The wider interior alone is not a convergence fix.
A targeted CPU-only saved-field diagnosis in job59953978 shows half-step bulk
H/M3.20227e-8/9.22185e-8, excluding observer truncation as the main cause.
Removing individual high modes5--8 does not substantially reduce the error;
in particular,7/8 have negligible effect. Filtered fields are diagnostic only
and are never accepted or exported as solutions. All evidence is preserved.

## Geometric interior extension trial — 18 October 2026

The published attenuated connection g*C is intentionally nonmetric-compatible
with the blended operator metric delta+g*(metric-delta). A new opt-in
HISPID_GEOMETRIC_TRUMPET_INTERIOR uses the Levi-Civita connection of that blended
metric for trumpet seeds only. This is a different interior extension; it does
not change the exterior equations, the flat plateau or physical seed fields.
The scalar operator uses its geometric Laplacian and the vector operator its
geometric div(L), restoring the latter's integration-by-parts energy identity
for compactly supported test fields. The coupled nonlinear system is not
claimed self-adjoint. Default behavior and the QI construction remain unchanged.

At three transition points the old metric-compatibility defect reaches0.560579;
the new normalized defect is1.82181e-16. Plateau/exterior cached operators agree
exactly. A separate-process API witness preserves QI residual/JVP bitwise;
centered trumpet JVP errors are1.34e-11 and8.10e-12. The small changed-interior
nonlinear control converges5Newton/57Krylov. These are implementation controls,
not physical binary acceptance.

The trial build enables the already-qualified batched modal transfer to reduce
factor-setup cost; other solver options stay fixed (GMRES64, forcing0.1,
nonlinear1e-12). One shared80GiB GPU allocation is requested for a small CPU/CUDA
field comparison followed by the192/224/240 study. The first level uses the
old192 checkpoint only as an explicitly migrated initial guess, with a fresh
solve; no old acceptance is inherited. Source/build suffix is
geometric-interior-20261018; output geometric-interior-sequence-v1. Driver and
local controls are in validation/trumpet/geometric-interior/. Physical outcome
is pending. The earlier wide-only failed sequence remains intact.

Job59954146 (session2285), one shared80GiB GPU for30minutes, builds and
qualifies the new image. CPU/CUDA changed-interior sampled fields agree to
2.30263e-13. The structural metric control also passes on the compute node.
Pinned producer8aba192865cd2b0b5be735bdbd6a3eef9f79a04a0f691ba56cb2cfab8b15e5db;
TwoPunctures dependency358c5e11907eca7037fe51a2c5e9f4a05a7ce5b126e9a21c512087b8e47f6331.
The192 solve is underway; no physical acceptance yet. A gated consumer script
is prepared for the240 solution, reusing the qualified CPU sampler and AthenaK
image only after fresh migration proof and three-grid acceptance.

The geometric-interior192 solve completed8Newton/50Krylov,82.9387solve seconds,
but fails physical momentum bounds: nearH/M RMS6.74833e-7/1.54108e-5 and
bulk2.57005e-8/1.11548e-5. The job stops before224/240. It is not selected
as a convergence remedy. Evidence remains in geometric-interior/evidence192.

A saved-polynomial interpolation diagnosis is now underway, without another
elliptic solve. The nodal differential operator and coefficient-based physical
sampler can acquire different roundoff tails when auxiliary low-mode nodal
values reach hundreds. The opt-in HISPID_BARYCENTRIC_SAMPLER uses tensor nodal
interpolation with constant-preserving difference sums and extended host
accumulators. It applies only to trumpet sampling; equations/JVPs and QI
sampling retain their previous paths. On the small control, residual/JVP
witnesses are bitwise identical and sampled-field difference1.57e-15.
An independent Cartesian-jet test uses a large constant plus highest-degree
radial/polar Chebyshev perturbations, checking vector correction values and
physical Kij at16/192 radial points. Locally, old/new Kij scaled errors are
7.76e-14/1.91e-15. The Perlmutter CPU control gives6.23e-16 for the new path.
This qualifies a sampler prototype, not a physical convergence result.
Job59954485 (session43649), a CPU-only15minute allocation, builds that sampler
and reobserves the saved wide224 and geometric192 bulk points. No other
compute allocation is active. The outcome must distinguish sampling roundoff
from actual discretization/solver error before another resolution is tried.

The saved-field barycentric recheck leaves the bulk momentum errors unchanged:
wide2249.22181e-8 versus9.22180e-8; geometric1921.11547928150e-5 versus
1.11547928157e-5. Sampler roundoff is not the dominant convergence failure.
Its prototype and controls are archived; it is removed from the production
source. The failed geometric-interior extension is likewise retained as a
replayable patch and removed from the production path. Neither was selected.

The frozen azimuthal probe59954609 (session43463) doubles16 to32 samples,
preserving the continuous field by orthonormal Fourier renormalization, with
no solve. In a meridional bulk band, native physical-equivalent residual RMS
at the old phases isO(1e-14); interleaved phase momentum RMS isO(1e-12),
maximum4.97e-11. These are far below the independent fully off-grid bulk
momentum errorsO(1e-7). An independent Cartesian FD witness agrees with
the interleaved native momentum components; Hamiltonian FD roundoff reaches
2.14e-9. The dominant remaining issue is meridional discretization/interpolation
between collocation nodes, not azimuthal resolution or the sampler.

Next candidate: dealiased Galerkin restriction of an overintegrated radial/polar
operator, retaining the original continuous CTT equations and interior
prescription. Preserve the original moderate observer profile for its new
physical sequence. Do not change physical tolerances or accept the failed
wide-only/nonmetric-compatible comparisons. All compute allocations from these
diagnostic trials have completed/released. No new Galerkin solve has launched.

## Overintegrated meridional candidate — 18 October 2026

The polynomial transfer test passes interpolation/round-trip/quadratic alias
errors2.60e-13/1.67e-13/6.22e-15. Projecting the weighted equations first
stalls the small linear control (relative0.388 after2400 iterations).
Projecting raw equations and applying N-grid row equilibration afterward
converges5Newton/3049Krylov,37.72CPU seconds. This is slower than the old small
control and is not a selected production improvement. The projected analytic
JVP has centered errors7.40e-12/5.46e-12; QI residual/JVP remain bitwise equal.
The prototype retains the original0.15/0.5 spherical interior and observer.
Quadrature is ceil(3N/2) independently in the two meridional directions, with
unchanged azimuthal resolution and N-grid tau conditions applied afterward.
The same continuous equations are evaluated; this changes the discretization.

The GPU path keeps tensor transfers, Q-grid fields/residual/JVP and the main
Krylov controller on device, with extended-precision host seed geometry.
The Q context never solves or builds factor banks. Kokkos's memory estimate
no longer double-counts a host PK_solve Krylov bank; host outer vectors and
both setup and resident Q-cache peaks remain bounded. No memory guard is
weakened or bypassed. The reference execution estimate is unchanged.
Job59955085 fails CMake generation from omitted test sources before numerical
work. Corrected job59955120 builds/runs one CPU/CUDA small control comparison,
then only the192×384×16 physical pilot. It uses the old192 original-formulation
checkpoint only as a migrated initial guess and must pass a fresh solve.
The candidate and all controls are in validation/trumpet/galerkin/. No physical
acceptance, speed improvement, finest-grid horizon result, or multigrid benefit
is claimed. The pending production study remains192/224/240 with one frozen
image and unchanged physical data/observer/tolerances if the pilot justifies it.

CUDA small control59955120 passes5Newton/3131Krylov,6.758solve seconds, with
CPU/CUDA field scaled difference5.25419e-13 (<1e-10). The default/OFF trumpet
residual/JVP witness is bitwise identical to the prior compensated image.
Native producerbff0cb239d8698f1be586c734671c946ec2b5bb088f7ac4d807ce2e9699e8d87.
The192 physical pilot is pending; no improvement claim follows from these controls.

The overintegrated192 pilot fails its first linear solve after2400 iterations
(444.507solve seconds); no Newton update is accepted. Setup211.720s,
total768.299s; job59955120 completes14:38 and releases. The original warm-start
field's physical errors remain unchanged; this is not a physical improvement.
The candidate is rejected and no224/240 or horizon run follows. Exact failed
result and raw observer are retained in galerkin/evidence192. Its option stays
OFF; the default/OFF exact equation witness remains valid.

A smaller preconditioner candidate addresses the retained grid-scale stalled
error without changing equations: two stationary residual corrections with
omega0.4 using the existing global modal FD inverse. The same generic helper
instantiates CPU and GPU operations, with three preallocated scratch vectors.
Default/OFF and QI retain the original undamped single inverse. Work counters
include both inverse calls and the extra Jacobian application. Native metadata
identifies the resolved preconditioner; trace mode prints recurrence progress
every100 iterations while true-residual confirmation remains mandatory.

The small CPU control converges5Newton/33outer Krylov (0.454s), with equations
and JVP bitwise identical to the old operator; converged fields differ1.68e-13
against1e-10. QI compatibility is bitwise exact and centered JVP errors are
1.34e-11/8.12e-12. Counts alone do not show a speedup: each preconditioner call
now costs two inverses and one extra JVP. No coarse-error bottleneck has been
demonstrated for multigrid. Job59955507 builds and checks one CUDA small control,
then runs one original-formulation240-grid pilot from the actual retained
failed state, capped at512Krylov per Newton. If useful, one frozen image will
be assessed at224/240/256, then its finest solution checked in AthenaK. The
corrected Kokkos memory accounting permits the planned256/restart64 level
within the same64GiB budget. No new performance/backend matrix is launched.
The attempted CPU-only overlapping build step on59955120 could not launch
(required ports in use); no compiler or numerical work ran there. The fresh
allocation supplies the first CUDA build for this candidate. Files and plan
are in validation/trumpet/modal-relaxation/.

The damped-modal240 pilot also fails:512outer steps, true relative0.77999236
against0.1, no accepted Newton update. Solve189.430s,total447.324s; job59955583
completed/released. Physical observer and original failure flag remain retained
in modal-relaxation/evidence240. No224/256 or horizon run is launched. The
option remains OFF; small-control success did not establish fine-grid benefit.

An opt-in small-host-matrix Arnoldi recorder now preserves unrotated Hessenberg
columns before Givens rotations. With the environment variable unset, it adds
no allocation or numerical operation. Existing independent Krylov controls
pass227 checks; diagnostic-off/on output is byte-identical. Independent least
squares on121 recorded small matrices matches recurrences to1.12e-17 scaled.
One128-step original-M retained240 replay is prepared to inspect the actual
preconditioned operator before another algorithm choice. The existing corrected
memory estimate permits restart128 at240 within64GiB. It is one bounded
spectral diagnostic, not a new physical or performance campaign. A changed
iterate is exported only as a diagnostic checkpoint if an update is accepted.
Files are in validation/trumpet/krylov-spectrum/. No physical gate is relaxed.

Arnoldi job59955993 completes one original-M128-step cycle, true relative
0.54673415, unchanged iterate. It reaches the previous long-run plateau within
one cycle, so further restart growth alone is not justified. The projected
operator norm is4.29e12 with singular-value spread about5.3e15 and large upper
couplings. Tail singular values and individual Ritz values are rounding-sensitive;
this does not establish exact full-operator eigenvalues. An untruncated QR
reconstruction matches the recurrence to7.9e-10 in relative-residual units;
ordinary SVD truncation drops directions and must not be mistaken for recorder
error. Raw data, analysis and scientific plot are in krylov-spectrum/.

The next candidate is left-preconditioned GMRES with explicit original-system
residual confirmation. A common backend-independent policy performs inner
M^-1 A cycles, updates the original step, and accepts only the original
unpreconditioned L2 target. The nonlinear1e-12 bound and all physical gates
remain unchanged. Local control passes5Newton/75Krylov,0.587s, residual/JVP
bitwise equal and converged sampled fields1.15e-13 from baseline. This is not
a fine-grid benefit claim. The descriptor records the side; native work
counts include the extra original-residual checks. Inner controller timings
refer to its composed operator and identity preconditioner, so whole native
solve time is the comparable performance measure. The bounded240 pilot runs
only after one CUDA small control. Default/OFF and other seed/system behavior
remain unchanged. Scripts and assumptions are in left-gmres/. Failed pilots
now optionally skip redundant physical observers while retaining their state.


## Regularity steering and QI control (18 October)

Left preconditioning59956254 failed at512 inner iterations: the explicitly
confirmed original linear residual increased to119652.944 times its starting
norm. No Newton step was applied. Failed checkpoints and small-control success
are retained separately; this is not a convergence improvement.

The user permits a different conformal-metric family if trumpet regularity
contributes. The fractional Schwarzschild exponent derived above applies to
the unmodified scalar sector; the implemented flat plateau removes that local
potential. It therefore cannot alone explain the plateau binary failure.
The nonanalytic compact transition, axis representation and modal conditioning
still merit separate consideration. A bounded QI Kerr control now keeps the
generic masses/spins/boosts/separation fixed, while scaling the interior window
to the QI throat. It uses the original modal inverse and existing CUDA image,
with a frozen160/192/224 sequence gated by the160 pilot's physical bounds.
This is a distinct initial-data formulation; no trumpet qualification carries
over. Near-hole observers follow each family's radius, so their norm ratio
is not a controlled slicing-only measurement. The bulk points and within-family
convergence can be compared. See regularity-family-control/ for the frozen plan.
No new multigrid or backend sweep is justified by the current evidence.


The narrow QI160 comparison59956570 fails the independent bulk momentum bound
(3.80871e-6), despite8Newton/46Krylov internal convergence. A wider QI interior
window, outer fraction0.85 with plateau0.15, reduces that error at the identical
exterior points to1.04547e-7. Its near H/M are3.13058e-7/9.14461e-8 and bulk
H6.30273e-9. One frozen160/192/224 study is running in59956729, with all prior
criteria unchanged; this candidate is not yet accepted. The wider window does
not require moving the QI observer shells. Actual binary horizon enclosure
remains mandatory. This is a distinct QI fallback, not a qualified trumpet fix.

An offline coefficient diagnostic using retained160 fields reveals slow m=0
meridional tails in both families. Degree8 polynomial endpoint windows leave
the trumpet wide-case radial vector tail almost unchanged when suppressing the
inner radial axis (ratio1.003), but reduce the selected degree80--143 band to
0.00406 when suppressing compactified infinity. The narrow QI case shows the
opposite localization. The inner radial boundary includes the whole axis
segment between the holes, not just the foci. These representation diagnostics
do not identify a logarithmic exponent or constitute a controlled slicing-only
comparison: the two retained interiors differ. The figure, source and input
hashes are in regularity-family-control/. Necessary low-m corner conditions
for a Cartesian-analytic correction are also recorded, not imposed in a solve.
For r=m modes they follow from u_m/ρ^m proportional to(1-t)^(m+1)P_m and
x=b(1+t)eta/(1-t): at(t=0,eta=side),
(d_t-2 side d_eta-(m+1))P_m=0. No extra tau experiment is launched while the
fixed wider-QI study is pending.


The fixed wider-QI160/192/224 study59956729 is complete. All three native
solves meet1e-12 and all declared absolute exterior bounds pass. Near-hole
H/M decrease, but bulk H/M increase to5.324e-8/1.699e-7 at224; the original
strict trend gate remains failed. Charge changes are<2.36e-9 scaled. The
combined66-point physical RMS decreases2.835e-7 ->1.406e-7 ->1.033e-7,
a supplemental summary that does not replace the failed regional criterion.

Finest224 sampler migration in59957080 agrees bit for bit at2454 points. The
first AthenaK attempt exposes a stale Python proof whitelist and stops before
the finder. AthenaK commit e54102e5 admits the explicitly bound v4/v5 bases
and rejects invalid maps/unknown versions; four regression controls pass.
Job59957101 reuses the passed migration without repeating native sampling,
and passes all horizon orders8/12/16 plus independent quadrature. Import
error is4.083e-16. Final MH=0.615673713/0.416086392,
Mirr=0.604595771/0.409618329, coordinate chi=0.370895278/0.345813441.
Final expansion RMS=2.041e-8/5.626e-9; buffered continuous enclosure margins
0.03205/0.01867 remain positive. Across all consecutive angular/quadrature pairs, maximum mass/spin changes
are2.005e-12/2.781e-11. These angular checks do not estimate the continuum input
error; the data and measurements remain diagnostic. No full trumpet or strict
QI convergence qualification is claimed, and no multigrid is added.

The PRD-style source and PDF include these results. Raw failed/successful
consumer attempts, migration witnesses, inputs, logs and surfaces are retained
in regularity-family-control/horizons/. A local retention audit verifies all
26 recorded horizon artifact references, both migration witness hashes and
all three solve-result/observer hashes. The broader convergence objective
remains open; this focused family/regularity follow-up and requested finest
horizon measurement are complete.


## Compactified-infinity and transition investigation (18 October)

All binary candidates in `validation/trumpet/infinity-regularity/` retain
`qi-moderate-wide`: masses .6/.4, separation6, generic seed spins/velocities,
exact seed scalar factor, same interior window, no far attenuation. Coordinate
maps are numerical discretizations, not changes to those free data. The
original three-grid absolute and strict near/bulk trend gates remain unchanged.

The host background now uses the existing algebraically normalized seed metric,
`expm1` attenuation complements and compensated trace evaluation. On an exact
constant transverse Schwarzschild conformal metric, the old quotient evaluation
introduced derivative error3.34e-6 on Perlmutter (.00694 locally) at tiny radius;
the stable evaluation gives
zero. The original-map binary prefix still fails bulk-H monotonicity, so this
fix alone is not a demonstrated cure for the physical trend.

A frozen-field doubled-azimuth check is far smaller than the observed bulk
error; raising phi is not justified. On the retained original-map224 field,
polynomial endpoint windows localize the scalar radial tail toward infinity
and the vector tail toward the inner boundary/transition. These are coefficient
localization diagnostics, not physical constraint norms or a measured
logarithmic exponent.

Three opt-in maps have distinct checkpoint identifiers. All are off by default:

* `HISPID_INFINITY_POWER=2` (v6) squares the compactified infinity distance.
  Its first160-grid linear solve failed; no finer solve was attempted.
* `HISPID_RADIAL_POLYNOMIAL=ON` (v7) uses `t=lambda*s+(1-lambda)*s^2`.
  This is unit-tested only, not qualified by a binary solve.
* `HISPID_DUAL_ENDPOINT_MAP=ON` (v8) integrates a positive polynomial density.
  With lambda=.05 and polar stretch3,
  `t=lambda*s+b*s^5+c*F(s)`, `b=(2-lambda)/5`, `c=1-lambda-b`,
  `F=s^15*(136*(1-s)^2+17*s*(1-s)+s^2)` (Beta(15,3) CDF).
  Endpoint slopes are .05 and2, preserving first-order infinity compactification
  while increasing density at both measured endpoints. The initial160 grid is
  cold; finer grids prolong the preceding same-map field. Acceptance is separate
  from the derivative and sampler checks.

The selected implementation uses host extended-precision geometry, CUDA GMRES,
16-thread ordered CPU physical sampling and OpenMP AthenaK horizon evaluation.
Parallel sampling warms the immutable spectral cache before starting threads;
all returned fields have a bitwise comparison to serial sampling. No multigrid
is adopted without evidence of a coarse-error bottleneck.

### Liu high-spin chart

`seed_family='liu'` implements Liu, Etienne and Shapiro, PRD80,121503 (2009),
[arXiv:1001.4077](https://arxiv.org/abs/1001.4077), Eqs.12/14:
`r_BL=r*(1+r_plus/(4*r))^2`. The isolated coordinate horizon is
`r_H=r_plus/4`, tending to `m/4` as chi tends to1. This is a wormhole chart,
not the existing R0=M trumpet. Exact extremality is still rejected. It reuses
the spin-axis rotation and full spacetime boost machinery; the radial metric
pullback and factored throat-limit curvature are evaluated analytically.

Independent QI-to-Liu tensor pullbacks agree to6.18e-15, including both bridge
sides, throat and axis. Finite-difference vacuum checks cover chi0,.39,.99 and
zero/generic moderate boosts. AthenaK reads the family and uses the matching
trial radius. An exact isolated, unboosted chi=.99 checkpoint passes sampler
migration and horizon orders8/12/16: mass errors below5.87e-13, spin errors
below1.37e-12, expansion RMS below9.84e-8. Horizon location has lower accuracy
than stationary mass/area integrals; no exact-extremal or Liu-binary validation
is implied. The fixed moderate binary study continues to use QI.

### Matched Bowen--York comparison

`validation/compare_fixed_by.py` matches collocation counts, bare masses,
separation, isolated Lorentz-transformed input spin and `P=m*Gamma*v`.
It uses exactly the stored HiSpID Cartesian observer points and FD steps.
Both physical Hamiltonian norms are `R+K^2-Kij*K^ij`, in the same units and
with equal-point RMS. Native weights, maps and free conformal metrics differ;
this is not a match of binary horizon properties or of internal residuals.
Failed stopping gates are retained and labeled, even when physical observers
are available. The capacity extension changes only checked limits/budget
plumbing: an existing supported BY solution and all exported derivatives remain
bitwise identical. The legacy allocation API retains its4GiB ceiling.

The historical default-QI sampler bitwise fixture predates stable background
arithmetic. The current original batch-seed exports remain bitwise identical
in516 values for QI/trumpet; numerical background/cache/sampler equivalence is
separately tested and may change roundoff. Physical free data and equations are
unchanged. This qualification does not weaken BY baseline preservation.


### Focus regularity candidate

The dual-endpoint map failed the first absolute bulk-M gate (2.143e-5),
although its nonlinear solve passed. It was not refined. A new candidate keeps
the original Mobius.2/polar2 map and QI-wide free data fixed and imposes
`(d_t-2*side*d_eta-(m+1))*P_m=0` at the two foci for modes0..4. Inside the
flat plateau, bounded Cartesian corrections solve a constant-coefficient
elliptic system, hence have removable singularities. This condition excludes
non-smooth low-mode corner behavior; it is separate from high-mode axis tau.
The native operator and directional JVP replace exactly two low-mode rows.
Normalized endpoint rows have coefficient L1 norm1; a rank-two boundary solve
makes the same row replacement in the approximate inverse. Moment reductions
and inverse corrections execute on the selected Kokkos backend.

A smooth Cartesian rational field and a symbolic polynomial give normalized
compatibility below9.24e-17 at96/144; an independent full dense inverse agrees
to3.65e-15. A small coupled CPU solve passes in9 Newton/32 Krylov steps. The
three-grid physical sequence is a separate obligation, declared in
`focus-plan.json`; no acceptance transfers from these algebra controls.

Focus tau solves reject non-flat interiors, nonpositive inner_min for an active
hole, and nonzero far_radius. The analytic far-reference correction W can be
non-smooth at boosted/trumpet punctures; homogeneous focus rows constrain P
alone and are qualified only when W=0. Sampling-only and raw operator contexts
remain available. The fixed QI-wide case already satisfies these restrictions.
GPU endpoint reductions use compensated FMA products and two-part partial-sum
joins, with an exact dyadic cancellation control. The CPU sampler/consumer
representation is unchanged; its frozen image is reused through an explicit
producer-to-consumer migration proof.


### Second-order focus compatibility and observer calibration

The C1-only fixed160/192/224 sequence meets all absolute constraints but fails
strict bulk trends. Sixth-order differences at4 times the same fixed FD4 step
confirm the rise, with bulkH4.38e-10,1.20e-9,2.24e-9 and bulkM3.89e-8,
5.02e-9,7.77e-9. The original FD4 verdict is retained. This higher-accuracy
observer was declared after coarse inspection but before inspecting finer
results. Rational-stencil, manufactured-curvature and exact Brill--Lindquist
controls pass; the exact-dataH floor is about8e-11. Saved-field FD4 replays
bind each observation; no additional elliptic solves are performed.

Matched-count vanilla-equation BY controls use identical physical observation
points and input bare charges, but distinct free metrics/horizon properties.
Their FD6 nearH8.03e-11,8.67e-11,8.09e-11 and bulkH9.94e-11,1.03e-10,1.11e-10
are verifier limited. All native1e-12 stopping gates failed at an original
reference residual plateau about1e-10; these remain diagnostic comparisons.
The supported baseline's33 arrays/derivative exports are bitwise preserved
by the explicit capacity-budget extension.

C1 solutions satisfy their normalized compatibility to1.7e-15 but have
second-order defects .00276/.00379/.00322. The next candidate uses, at each
focus for modes0..4, C1 plus
`P_tt -4P_etaeta -2(m+1)P_t -4*side*P_eta +m(m+1)P=0`.
This follows independently from a smooth Cartesian field and
`Q=(1-t)^(m+1)P`, `Q_tt=4Q_etaeta+4*side*Q_eta`. Four rows replace the old
two; the exact rank-four inverse, row scaling, endpoint moments and correction
are shared across CPU/GPU. Checkpoint IDv10 separates this row system fromv9.
The physical formulation, original compactification and acceptance gates stay
fixed. No analyticity at infinity or multigrid benefit is assumed.

The first GPU preflight caught missing unit scaling on the two new rows before
any large solve ran. The corrected small coupled solve takes9 Newton/35 Krylov
steps, versus1164 Krylov in the rejected implementation. Smooth-field
compatibility and independent dense inverse errors are2.11e-16/3.24e-15 locally.
The frozen pure CPU sampler built before this solver-only scaling correction
is intentionally preserved: sampler equivalence will be established separately
by actual producer/consumer migration, not inferred from this algebra test.


Both focus orders are now selectable without duplicating the operator or
inverse: `HISPID_FOCUS_ORDER=1` is the default when focus tau is enabled;
`=2` selects the experimental v10 row system. `HISPID_FOCUS_TAU` itself
remains OFF by default. The C2 coarse160 run has bulkM1.053812e-6 and fails
the unchanged1e-6 bound. Its separate192 refinement gate is declared before
inspection; no failed result is hidden by a later choice of grid range.
The existing fixed C1 three-grid sequence and its FD6 failure remain the
primary record until another candidate meets all gates.


### C2 refinement rejected; retained C1 result

The C2 refinement check worsens bulkH/M to6.29684e-9/1.77824e-6 at192,
versus1.95912e-9/1.05381e-6 at160. Both internal solves pass (5 Newton,
35/33 Krylov). The unchanged physical momentum bound fails and the
predeclared improvement gate rejects refinement. No224/256 C2 trial is run.
This stronger tau condition is not adopted as a convergence remedy.

The retained first-order160/192/224 sequence improves finest bulkM by21.9
to7.77e-9 (FD6), and meets every absolute exterior bound, but strict bulk
H/M trends remain failed. Physical metric/Kij Cauchy changes decrease and
charge changes are below4.08e-10; neither replaces the failed constraint gate.
Direct current-order-one row evaluation on the saved states remains below
1.60e-15 at all grids; old/new normalized weights differ by at most3.47e-18.
The shared selectable implementation is therefore compatible with the
retained boundary states, without claiming a new physical solve.

No further backend sweep or multigrid experiment is adopted in this accuracy
study. The selected C1 grids take28/27/238 total Krylov steps: finest-grid
cost growth is real, but a coarse-error multigrid remedy has not been
demonstrated. All three collocation solves already reach tolerance while
the independent off-grid physical trend fails.
The remaining strict bulk convergence objective is explicitly incomplete.


### Finest retained C1 horizon completed

CPU allocation59965595 checks the actual224/448/16 checkpoint with the
frozen CPU sampler and AthenaK executable. Producer/consumer fields and
dgamma are bitwise identical on2454 points. All8/12/16 angular orders and
independent32->48 quadrature pass, with native geometry callbacks and zero
evolution. Final MH=.615673712968/.416086392223, Mirr=.604595770588/
.409618328822, coordinatechi=.370895278083/.345813441147, expansionRMS
2.040501e-8/5.625873e-9. Maximum mass/spin-vector refinement changes are
3.50e-13/4.04e-12. Modified-interior enclosure margins after the observed
shape buffer are .03204669/.01867393. The strict bulk constraint gate
remains failed, so stronger binary validation is still false.

`final-integrity.json` verifies all retained horizon artifacts, both sampler
witnesses and all three physical observers against their recorded hashes.
No compute allocation remains active. AthenaK commitde640365 is pushed.


### Lower-resolution study and direct single-hole grids (18 October)

The unchanged QI-wide/C1 formulation is now solved at8/16/16,12/24/16,
16/32/16,32/64/16,48/96/16 and64/128/16. The BY comparison matches bare
charges, collocation counts and the same66 exterior Cartesian observers;
it does not match the conformal free metric. FD6 uses the same previously
calibrated steps. All low-grid collocation solves converge. Near and bulk
sample H/M RMS decrease for QI across this sequence. For32/48/64, nearH is
8.9273e-5/4.2249e-5/1.6043e-5 and bulkM4.4965e-5/1.3166e-5/6.3237e-6.
These grids are below physical acceptance; they demonstrate the truncation
regime and do not reverse the retained high-resolution bulk-trend failure.

BY nearH decreases from.01210/.002801/.0004875 at8/12/16 to3.8039e-7 at32,
then reaches the observer floor:9.8858e-11/8.4016e-11 at48/64. Its64-grid
metric differs from retained224 data by9.09e-15 in max |delta|/(1+|reference|).
Vanilla line/BiCGStab and modal CUDA GMRES agree to5.37e-16 on the moderate32
grid. QI cold and prolonged64 solves agree to1.78e-12 or better in every
physical field. Independent solution-reference comparisons, including the
original compactified6^3-point U norm, use saved fields without new solves.

Paper convergence benchmarks use single-hole S/m^2=.2 or P/m=.2 and binary
m+=m-=b=1, P+=-P-=.2 along y. The single-hole BY spin relative-U error falls
from2.59e-7 atN10 to2.31e-12 at16 and then roundoff; the boosted case has a
fitted order3.976 overN10/16/24/32. Its n=70 reference and spin n=50 reference
use the paper's parameters. Legendre angular collocation replaces theta
Chebyshev, so these reproduce the physical equations and convergence test,
not every discretization detail. The binary uses16/24/32/48/72 with even
Fourier counts;72/72/36 replaces the paper's70/70/35. BY solution differences
fall rapidly and the infinity error becomes dominant. No accurate fitted
sixth-order asymptote is claimed from these few points. The72-grid native
1e-12 stopping gate fails at7.95e-11, while its sampled fields agree with48
to near machine precision. That failure is retained, not relabeled converged.
The Kerr binary counterpart has separate, fully specified QI free data and
shows decreasing solution differences and constraints. The old BAM fixed-mesh
multigrid comparison and test-mass binding-energy study are not reproduced.

`python/single_hole.py` adds a direct two-dimensional CPU path with
r=(m/2)A/(1-A), Chebyshev radial roots and Gauss--Legendre polar nodes.
The BY nonlinear solve uses a Legendre-block inverse of the flat operator;
the Kerr path stores regularized spherical tensors from native exact seeds.
Both support axial spin and momentum/boost, and no inactive focus. Existing
native binary layouts/defaults are unchanged. Direct BY data agree with the
native binary-map single-hole control to1.10e-15 (spin) and1.10e-12 (boost)
at the stated points. Independent nonzero initial guesses agree to4.17e-17
in the direct BY unknowns. This is not a perturbed-Kerr coupled elliptic solve.

Liu chi=.9999 is represented at64/32,96/32,192/32. Scaled physical metric
and Kij differences from exact data decrease to1.60e-13/2.98e-12 on14 fixed
points, including both sides of the horizon. Independent H/M RMS decrease
to1.31e-8/6.87e-10, matching the exact-seed FD6 floors on those same points.
Known-sphere horizon integration gives M=1.0000000000000007 and
chi=.9999000000000002; expansion RMS is below7e-13. Independent radial tensor
pullback agrees to7.90e-14. The FD4 vacuum check needs.002/.001/.0005 steps
at this spin: the coarse .001 result exceeds its unchanged1e-7 bound, and
.0005 passes with the expected fourth-order decline. This is analytic-seed
and representation convergence, not an AthenaK chi=.9999 finder or evolution.
The Liu chart remains a wormhole with rH=.25353544551648M here.

All binary work used one shared interactive A100, with16 host threads; the
small spherical problems use CPU NumPy/SciPy. Selected results, source/image
hashes and plots are under `validation/trumpet/lower-resolution/`; larger
unknown arrays remain on Perlmutter in `lower-resolution-20261018-v1`.
Scripts record the actual stages and library identities. A saved-norm replay
caught an archived shared-library SONAME substitution when BY was loaded first;
the replay now loads and verifies HiSpID before the explicit BY image, as the
successful solve driver already did. This failure involved no physical solve.

A single spectral shell is plausible for linear black-hole perturbations,
but this elliptic grid is not yet a time-evolution discretization. For an
exterior evolution, a horizon-penetrating slicing/excision boundary and a
finite outer boundary or hyperboloidal compactification are separate choices.
Liu wormhole slices and compactification toward spatial infinity cannot simply
be assumed to provide those conditions. General perturbations also need
angular modes, physical free-data specification and new constraint solves.
No work on Lazarus/Teukolsky evolution is performed here.

## Physical metric/time-derivative free-data interface (2026-10-18)

`HiSpID_create_with_free_data` reuses the coupled CTT solver with a per-context,
reentrant callback supplying physical lab Cartesian delta-metric and fixed-
coordinate delta-metric time-derivative jets. The reference stationary Kerr
lapse, shift and conformal factor define the ADM/CTT split. Metric Hessians
and gdot gradients are required; no finite-difference cache construction is
introduced. The callback remains alive until context destruction.

`HiSpID_create_with_wave22` provides a native pole-regular electric l=m=2
packet matching SpEC commit 6b55d93e51076dc75a006837e808f1182549bf64. A stationary
smooth radial cutoff protects the throat. The supported reference is exactly
one unboosted active QI/Liu hole; the existing full 3D two-focus map and four
unknowns are retained. The disabled hole still defines a distinct map focus.
All old ABI layouts/constructors retain their meaning and arithmetic.

```python
from hispid import Backend, Hole, Wave22
b = Backend('/absolute/path/to/libHiSpID.so')
c = b.config()
c.hole[0] = Hole(1., center=(0,0,0), spin=(0,0,.9))
c.hole[1] = Hole(0., center=(0,0,-12))
c.seed_family = 'liu'
c.n[:] = (96,192,16)
c.memory_limit_mib = 32768
c.tolerance = 1e-12
c.max_newton, c.max_krylov, c.krylov_restart = 12, 1200, 100
c.far_radius = 0
c.omega[:] = (0,0)
c.inner_min[:] = (.1,0)
c.inner_max[:] = (.2,0)
with b.create(c, execution='kokkos', wave=Wave22(amplitude=.001)) as s:
    diagnostics = s.solve(krylov='gmres')
    # Save this context's actual parameters with its nodal unknowns:
    if diagnostics['converged']:
        s.write_checkpoint('/absolute/path/to/packet.hispid')
    # A nonconverged diagnostic remains sampleable for diagnosis.
```

Host long-double geometry remains the qualified setup stage. Its cache is
uploaded to the unchanged production CUDA solve. Arbitrary Python callbacks
acquire the GIL during setup/sampling; heavy production providers should use
the native callback API or the built-in packet.

Packet checkpoints use explicit version3, including family and nine wave
parameters. The corresponding AthenaK reader reconstructs the same free data.
Generic callbacks cannot be silently exported through the old checkpoint:
applications must retain their own provider representation.

The trial time derivative is not fixed physical gdot after projection. Rebuild
gdot from the solved metric/K, lapse and shift to inspect its radiative modes.
This feature supplies vacuum gravitational-wave initial data, not scalar-field
data or a superradiant evolution/instability calculation.

Cold starts and three-resolution results are under `validation/trumpet/wave22/`.
The inactive focus at z=-12 improves packet resolution at fixed physical free
data. A deliberately incompatible constant-P warm start fails; it violates
the homogeneous focus rows. Use zero corrections, compatible prolongation,
or an admissible correction guess. This failed diagnostic is retained.


## Near-threshold perturbed Liu chi=0.9999 (18 October2026)

The original prograde/ingoing electric l=m=2 packet passes the finite-observer
checks on the production CUDA radial/polar sequence96x192x16,128x256x16,160x320x16. Its mass-one seed
spin is.9999, amplitude.001, carrier omega=.99(2 Omega_H)=.9760972454792889,
cycles4, r0=12 and cutoff2--4. The inactive focus is fixed at z=-32.
The finite packet spans the superradiant threshold in frequency.

Independent54-point FD6 Hamiltonian RMS decreases5.0492e-8->2.7922e-8->1.6498e-8
and momentum7.2817e-7->4.0719e-7->2.6097e-7. Physical solution differences
decrease; a compatible1.05-scaled start returns sampled metric/K within
8.7e-12. The observed rate is about second order, with fixed Nphi16; neither
global nor exponential convergence is claimed. A focused regression fixes
the time-derivative observer's old hard-coded chi=.9 gauge assumption.
The shared Liu configuration/gauge lives in validation/kerr_wave_case.py,
independent of the oracle test. The study CLI records actual spin, carrier
threshold ratio, packet parameters, radial derivatives and source hashes;
nonconverged checkpoint exports retain their diagnostic label.

The finest portable checkpoint passes bitwise9-field CPU migration at1602
points and AthenaK's lmax8/12/16 plus32->48 quadrature checks, unchanged
expansion RMS1e-7 and inner-ball enclosure. Final MH=1.000000000118978,
Mirr=.7120891092880841, chi_coord=.999899999762024, expansion RMS9.8172e-8.
These are actual direct-native-geometry finder searches with zero evolution
steps. Coordinate spin is not AKV spin, and horizon mass is not total ADM
mass. The final asymptotic metric-characteristic outgoing/incoming amplitude
ratio is.00269 at r=24.8741: predominantly ingoing, not a curvature-flux or
amplification measurement. Source, checkpoint and loaded-image bindings,
raw traces and one compact figure are in validation/trumpet/wave9999/ and
the existing PRD report. The earlier .9 campaign and failures remain intact.


## Physical radial-resolution audit (19 October 2026)

The original chi=.9999 packet grid is under-resolved. Its fixed-point and
horizon results above remain valid at those observers, but do not qualify
the packet waveform. Physical radius changes both prolate coordinates;
increasing the radial count while reducing the polar count is insufficient.
`validation/wave22_physical_radial_audit.py` interpolates the complex carrier
on the actual two-dimensional grid and differentiates along nine physical
radial rays, including one-sided probes around the inactive focus. Errors
are normalized to the wave peak, and phase is checked above 10% peak.

At the old160x320x16 grid, value/derivative errors are29.5%/42.9%. With the
same physical free data, inactive focus z=-16 and radial stretch1, the
576x128x16 grid passes0.1% value/derivative and.001rad phase gates.
The new map needs its own matching producer/CPU reader: the continuous
basis identifier ends in `_r1_k2`; the old `_r0.2_k2` reader cannot read it.
The solved-wave profiles and strict native stopping gate are separate
requirements. Nonconverged diagnostic checkpoints are rejected by the
profile qualification tool. The original evidence remains unchanged.

The restart200/4000-iteration production continuation also fails: native
scaled maximum9.57e-11 versus1e-12; the second inner solve stops at true
relative L2.642 after4000 iterations. H RMS7.97e-10, M RMS2.12e-8 at54
observers; M is worse than the earlier diagnostic3.90e-9. The experimental
spectral radial inverse is not promoted. No revised n480/n576 solve or
new AthenaK horizon passes are claimed; those are carrier-only audit grids.
The original-data dense wave profile check measures both conjugate m=+/-2
modes, with all differences normalized per angular ring to the wave peak.

Fourier localization of the final diagnostic puts99.56% of squared scaled
residual norm in m0--3. Their maxima are at r=.1745 (interior transition)
and r=.294--.301 (just outside the seed horizon). Parseval disagreement is
5.8e-15. This localizes native equations, not physical constraints, and
does not establish a cause or justify any stopping-norm relaxation.

Dense profiles of the original saved96/128/160 solutions fail the0.1%
wave gate. Finest-pair peak-normalized time-derivative difference=.1047114,
ingoing=.0523327, metric=1.8718e-6; differences decrease but remain too
large. At485 packet FD6 points (97 radii, five rays including both poles),
finest H/M maxima6.36e-8/6.86e-7 fail1e-8. Raw seed calibration maxima
5.44e-10/3.13e-16 are retained; no floor is subtracted. Full checkpoint
headers confirm the same formulation. Updated report and one compact
figure distinguish old solved-data failures from proposed carrier-only
grids. Job59975943 is released; no revised horizon search was run.


## Curved inverse and exact polar storage (19 October 2026)

`HISPID_CURVED_MODAL_PRECONDITIONER` is opt-in. It retains the cached curved
scalar Fourier symbol and 10/9 times its vector principal symbol, including the
strong mixed radial/polar principal term. The nine-point stencil stores
three radial bands between polar blocks; geometric row equilibration moves
the reciprocal factor into the RHS scale. Residuals, Jacobians, free data,
physical boundary equations, and stopping criteria are unchanged.

`HISPID_STREAM_POLAR_PRECONDITIONER` requires the existing exact polar
correction. It batches the Schur construction and applies the same correction
through a second base inverse, avoiding stored dense response banks.
Independent stored/implicit full-matrix witnesses agree below 7e-15; actual
CUDA fields/JVP/inverse controls pass. At 576x192x16 the removed responses
alone are 7.59 GiB per host/device bank; this is a storage count, not an
end-to-end peak-memory measurement. The bounded 32x64x16 control preserves
residual/JVP bitwise and returns the same sampled geometry within 8.9e-11,
with 52 versus 111 Krylov iterations. This grid does not resolve the packet.

At 384x128x16, low-mode FD factor witnesses are accurate below 9e-13 and
dominant CPU/GPU inverse actions agree around 1e-14 relative. Nonetheless
the equilibrated mixed/exact-polar production continuation fails after 2000
iterations at true relative residual 0.6432 (target 0.1), 225.43 s, accepting no
Newton update. It is not a convergence fix. The native maximum remains
2.05e-10 versus the unchanged 1e-12 gate. Evidence, including failures, is
in `validation/trumpet/wave9999-curved/`. Coupled multigrid trials did not
improve the production solve and were retired; their source and traces are
preserved there. A two-level physical trial failed its 0.1 linear target
at 0.870 after 600 steps. Subsequent transfer controls passed, but the
manufactured solve missed its unchanged target; a deeper hierarchy was
cancelled for cost without a production run.

The separately opt-in `HISPID_CONSTRAINT_BLOCK_PRECONDITIONER` solves the
full vector block, forms the exact Hamiltonian/vector coupling, then solves
the scalar block. The original coupled Jacobian and outer/native tolerances
are unchanged. Its first 64-step/restart-12 inner version passed the coupled
manufactured witness in 9 outer steps (439 vector and 41 scalar inner steps),
but its physical continuation was cancelled after 100 outer steps at 0.858
without a Newton update. This is not a completed 300-step failure.
The next candidate packs active block components, skips only exactly zero
modal RHS tasks, and uses restart 96/maximum 256 inner steps with the same
0.1 block target. Setup and solve memory peaks are bounded separately;
the device free-memory guard remains. Qualification is pending. No revised
three-grid or AthenaK horizon pass is implied by these algebraic controls.

The public radial extent cap is now 1024, polar 512, azimuthal 256, with
explicit aggregate memory budgets through 131072 MiB. The C configuration
layout is unchanged, and allocation/device-free-memory guards remain.
A focused actual AthenaK C++ reader/CPU sampler control passes 576x8x4
and rejects radial 1025 and an over-limit budget. A matching native header
and reader image are required; this is not a refined horizon check.


## Retained search vectors (19 October 2026, candidate V16)

The fixed-M single-basis recombination experiment passed small algebraic
controls but developed a large physical true-residual gap (first vector
block of the next Newton step: 2.335e4 after 256 iterations). It was cancelled
before its 40-outer-step pilot budget; its source is archived, not retained
as an additional public Krylov method. No solution is promoted from it.

The new CUDA-only storage policy retains each actual preconditioned search
vector in lossless HostSpace doubles. One device staging vector restores
the original FGMRES update, including variable approximate inverses. All
Arnoldi, A/M and true-residual calculations stay on device. CPU storage and
the default equation-independent solve remain unchanged. The 256-step inner
limit/restart and 0.1 target remain fixed; the outer/native targets remain
0.1/1e-12. Aggregate allocation bounds include the retained host bank.
Actual CUDA retained-Z controls pass bit-for-bit for the same variable-M
GMRES/LGMRES cases (restarts 3 and 19). At restart 19, device Krylov
workspace is 3648 versus 6384 bytes, with 2888 host bytes; this is moved
storage plus staging, not a whole-process RAM reduction. Core C controls
227 and CUDA controller/coordinate controls 286329 pass. The unchanged
coupled witness passes in 7 outer steps, original true relative residual
2.2918980572e-8 and correction-jet relative error 9.3545147046e-7. The
three-resolution physical study remains pending. Producer image SHA256:
fa97396c525bc8663d2b65d0fa80320a2aa7f68c5dcf8627bcd28fbe32935e54.
These algebraic controls do not establish the physical convergence fix.


The completed V16 outer-restart-12 physical pilot accepts one additional
Newton update (four from zero), then misses the original 0.1 linear target
at true relative 0.1011709103 after 40 outer steps / 10240 vector-inner /
306 scalar-inner steps. Native maximum 1.8934877815e-9 still misses 1e-12.
Native solve 1218.92265 s, excluding geometry 61.52597 s. No physical pass
is implied. A detached continuation is prepared with the SAME producer
image and existing LGMRES, outer restart 24; no compilation or controller
tests are repeated. Original tolerances/formulation are fixed. The finest
allocation bound including LGMRES headroom is 38.955 GiB guarded GPU /
91.544 GiB aggregate, below the 96 GiB context budget. The completed pilot
checkpoint is used only as an initial guess. SSH was lost before collecting
the final checkpoint receipt and new allocation state; neither is assumed.


### Recovered outer-24 continuation (20 October 2026)
The detached V16 LGMRES(24) run, job 60009100, hit its two-hour allocation limit before its first new Newton update. Its last logged recurrence residual was 0.10107067270049126 at outer iteration 140; there is no terminal true-residual result. This is not a completed 300-step failure. No revised resolution or horizon gate passed. The saved four-update diagnostic checkpoint has SHA256 a66e1d57ebda086eb851fbe14c1c0a11c898827f5ecbb1be88967c2a850f4d18 and independently sampled 54-point H/M RMS 2.5239566021653973e-9 / 2.387272605330194e-8. The full receipt and compressed interrupted-run log are retained in wave9999-curved/interior-focus.


### Inner recycling and recoverable continuation (20 October 2026, V17)
The frozen four-update physical RHS is dominated by momentum residual (99.9992071% of squared L2), not high-mode polar tau (6.18e-11 fraction). CPU/GPU M actions agree to 3.30e-13. Existing low-mode focus rows are separately reported in V17; this does not make the smooth manufactured witness a physical qualification.

The opt-in constraint-block inverse now uses the shared LGMRES controller, restart 256 / limit 1024 / three retained corrections, with the same original true 0.1 target. No equations, map, free data, public ABI, outer controller, default BY behavior, or acceptance thresholds change. Actual physical pilot: 490 vector + 3 scalar steps, one outer step, original true relative 0.09906689597, full accepted Newton update, 88.37384 s excluding geometry. Native maximum 2.26938e-10 still fails 1e-12; the diagnostic checkpoint now contains five updates from zero. Physical H/M RMS 2.51541e-9 / 2.64357e-8, with M slightly worse than before. Three-grid/horizon acceptance remains unproven. Producer SHA256 1eec375773a0fd8f287320b08972364e075441c35fb4b1c1623d48ceeff1a4b0. Local OpenMP6 and CUDA coupled witnesses retain seven steps / true 2.2919e-8 / jet 9.3545e-7. Finest outer-GMRES24 conservative allocation: 38.7044 GiB guarded GPU / 91.4120 GiB aggregate, including extra augmented directions and lossless host Z.

The optional Python checkpoint workflow calls the native solver for one Newton attempt at a time, with unchanged forcing and gates. Only native Newton-limit status permits continuation; linear/line-search failures stop. Each completed checkpoint and each JSON report is individually replaced atomically (they are not a joint transaction, so recovery verifies the actual checkpoint hash). A converged coarse solution can be spectrally prolonged only as an initial guess. Work/history across calls are aggregated. Native factor caches are rebuilt per call. Three targeted recovery tests pass. A fresh single-GPU shared-interactive job 60056582 runs the fixed 384/480/576 sequence and gated near-hole, packet, portable-reader, and AthenaK checks; launch is not evidence of completion.


### Frozen next-step diagnosis and coupled candidate (20 October 2026)
The V17 full continuation was cancelled during the second outer application, with no further accepted update; first vector true residual 0.99438137 after 1024 iterations. This is not a completed 80-step result. The saved five-update checkpoint remains diagnostic.

The isolated V18 full-axial-diagonal candidate completed its six-outer-step pilot on the same checkpoint: original full true relative residual 0.94842698568, 6144 vector + 6144 scalar inner iterations, 1357.84703 native seconds. No accepted update; the unchanged native maximum is 2.26938e-10. The source is reverted and archived separately. Job 60057063 MaxRSS 10733796 KiB. The full receipt and log are retained.

V19 is a default-off coupled meridional approximate inverse: axial Cartesian m=0 plus radial Cartesian m=1, preserving the complementary and toroidal components. Its coefficients use the fixed curved vector operator; GPU gather/scatter surround a cached host factor. Seven-band interleaved LU and component-specific focus constraints share the existing modal machinery. Equations, free data, maps, stopping norms, and public ABI are unchanged. Local cross-action, independent inverse, bridge-complement, and existing coupled-system checks pass. GPU/physical qualification remains pending. Perlmutter rejected the initial 64 GiB host request under its one-GPU/32-core shared-queue rule; job 60058420 uses 56 GiB. The aggregate bound is conservative and is not a measured host peak. The fixed three-grid and AthenaK pipeline remains gated on actual native acceptance.


V19 CUDA controls passed: cross-action 4.29519e-13, seven-band inverse 8.88e-16, existing full-system witness 7 outer / true 1.3193931858e-8 / jet 2.5715666133e-8, bridge 5.44e-17. Producer 218226736225f213f199d9dd592a7229133c8cda151b5b575cf76f0f5e0e4ac4. Physical first vector block still ends at true 0.99286642 after1024 steps (scalar89/true0.0929590). Job60058420 was deliberately stopped after this first completed block application, without an accepted update or terminal outer residual. This is not a completed six-outer-step failure. Candidate source is archived and reverted; production code remains V17. A one-outer-step replay of the exact V17 image, with existing raw-Hessenberg export and no arithmetic changes, is prepared to test whether discrete weak modes justify deflation. No native/three-grid/AthenaK gate has newly passed.


The V17 one-outer-step replay completed as job60059283: true0.99317869884,1024vector/129scalar inner,158.35325s; no update. First unaugmented256-column H norm2.41314e21, columnnormspread3.17216e14;161computed singularvalues beloweps*norm. DoubleQR cannot reproduce the native recurrence; tailSVD/leftnull estimates are unreliable. Decimal35/65digit Givens and backsubstitution agree at0.9926147260210148780680 on the exact recorded roundedH, so higherprecision projectedLS alone is insufficient. No deflation implementation is justified from this projection. RawH, scripts, receipt and log are retained. V20 adds only an opt-in no-Krylov modal-response audit, leaving production arithmetic unchanged. It reproduces the normalized first vector-inner RHS and base inverse, splits unknown modes0..4/ge5, applies the exact momentum block, separates PDE/focus/tau actions, measures physical collocation field/gradient norms, and checks superposition. This is diagnostic instrumentation, not an accepted solver fix.


V20 frozen modal-response audit completed, job60059703, producer8fc24511e0e7a020f3df126511b0901ea37ce58780a1b9541e703a3f3c0c68d4. Original vector RHS L2=9.5271679411e-9, normalized to1 for the diagnostic. Actions for groupsm0,m1,m2,m3,m4,ge5 are1552.30,597.56,26488.44,2564746.05,143629180.17,14.90. Amplification is in ordinary PDE rows; largest m4 response at r1.3732, not a horizon constraint measurement. Superposition relative3.88954e-16 and group RHS Parseval sum1. No solve or accepted update; last_gmres_relative1 is an explicit diagnostic placeholder. Field/gradient norms are unweighted Cartesian collocation sums. Result/log/recipe and structured summary retained. This supports targeting the approximate-inverse response; no numerical remedy or physical convergence is claimed.


V21 spectral-polar inverse completed diagnostic job60060245. The generic block solver was extended to dense inter-block coupling and transposed focus rows, with the fine A/formulation/basis/gates unchanged. CPU/CUDA polynomial, independent full-matrix inverse, np4/np16 bridge/Nyquist, and existing coupled controls passed. CUDA witness7outer,true2.1201016495e-8,field1.2231821330e-7. Producer a456017be0c401b59348462a4ba7ae0fa91b9c37e5b19f14f1ecd2b44139b6a2. Actual frozen combined action5.8493890579e7, ratio0.4071880858 vsV20; superposition3.90687e-16. Failed predeclared0.01 action-ratio diagnostic gate: no Krylov pilot/accepted update. Candidate archived and reverted. This is not a completed physical-solve failure; it shows exact polar differentiation alone does not sufficiently reduce amplification. No finest-grid allocation was attempted: extra bounds4.20GiBdevice/8.24GiBaggregate exceeded the previous resource plan. No memory guard or physical tolerance was relaxed.


V22/V23 frozen diagnostics localize the base-inverse amplification at compactified infinity. Actual m4 P correction has >0.9999999999 squared norm in final4 radial rows; maximum r3.8247e6/rho45564 (not physically near axis). Offline inverse-load reconstruction uses ordinary long-double Fourier projection, not exact native compensated summation. V23 exact m4 fine action before focus=143597441.1421, focus increment=39457.35936, combined=143629180.1652. No Krylov solve or update. Producer30abd0019c3c48625f4a4699d2dc1a24eb72a4b41de0c487aaa15f8d8c0d4279; original checkpoint unchanged. Raw arrays stay in local tmp/remote scratch with committed hashes, compact analyses and recipes.


### Rejected outer inverse and tail-fit controls (V24/V25)
The uncapped inverse-only damping was cancelled locally after 3 min 21 s without a terminal residual; this is a cost cancellation, not a completed accuracy failure. The capped follow-up, GPU job 60061511 / producer 9280d93516d2f4feefb933f556973ef4adf8883bdf21a12d04d8fc19688a3456, fails the unchanged coupled manufactured field gate: 16 outer steps, true residual 7.1313394703e-8 < 1e-7, field/first-derivative relative error 0.02886591537 > 1e-4, 93.07579 native seconds. The printed solver success is not a field-recovery pass. Neither candidate reached the physical pilot; both are archived and reverted. No formulation or acceptance threshold changed.

The retained five-update checkpoint has large outer modal-P components, but P amplitude alone is not a physical error norm because the physical basis vanishes at infinity. One fixed quadratic-in-inverse-radius tail fit over confocal radii 256--1024 fails alternating held-out samples for m=2/3/4 (relative P errors 1.55/1.95/2.00). No edited checkpoint or fitted initial guess was produced. Half-turn decomposition is also read-only; it motivates a cold-start symmetry-sector experiment, not projection of an accepted solution. Compact evidence and scripts are under wave9999-curved/interior-focus. Correction to the preliminary rationale: 24.8741 is the packet midpoint, not its end; the nominal outer window transition is 37.7482 with exponential tails, so the rejected damping onset at 32 was inside that window.


### Optional half-turn unknown sector (V26)
`HISPID_HALF_TURN_SOLVE` is default off and requires the modal axis-tau basis. For a single unboosted Liu hole with spin and grid aligned with lab z, the native real l=m=2 TT wave is invariant under a half-turn about that axis. Scalar and axial-vector modal unknowns retain even m; transverse Cartesian-vector unknowns retain odd m. Arbitrary callbacks and unsupported seed/boost/grid combinations are rejected. Only initial unknowns and inverse outputs are restricted; all original RHS, Jacobian output, nonlinear residual, line search, linear true residual and acceptance norms remain full. All original factor banks remain present. The axisymmetric equation projection cannot be enabled together with this option. A nonzero forbidden residual can therefore cause an honest failure rather than being hidden.

CUDA V26 control, job60061919 / producer b9517e9895f5f0a495e1656d73f2572a0806a6b4742705fbc5b983a715672bb6: TT metric/gdot jet covariance error0; unprojected zero-source parity1.81253e-16; full Jacobian parity6.13558e-17. The actual-wave-geometry 128x96x4 manufactured correction passes original true1e-7 / field-first-derivative1e-4 gates in7outer, with errors1.52721e-8 / 9.89321e-8 and0.98052native seconds. This is an algebraic control, not radial/azimuthal wave resolution. The later local guard/test forbidding combined half-turn and axisymmetric modes is absent from this immutable GPU image; the tested high-spin wave already independently forbids axisymmetric mode. Independent read-only review found no remaining blocker for the guarded candidate.

The first physical launch stopped before solving because an older Python driver did not support per-Newton checkpoints. Corrected job60062182 uses the same hash-verified library and study-wave-curved-v17 driver, starting from zero on384x128x16. Eight Newton attempts / six outer iterations bound this pilot; native1e-12 and original linear relative0.1 are unchanged. It is not a completed physical convergence result. Failed-candidate archives and controls do not qualify the pending three-resolution or AthenaK horizon study.
