# Trumpet implementation and acceptance ledger

Goal started 2026-10-05. Native base4ce5974; AthenaK consumer base9495fdb9.
Dedicated branches: codex/hispid-trumpet and codex/hispid-trumpet-pgen.
The prior worktrees and raw results are unchanged. No previous benchmark
campaign is resumed. Current status: seed-family API, context/cache/sampler
dispatch and versioned checkpoint integration implemented; no trumpet binary
or production-backend qualification yet.

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
There is no claim to computed K Hessians. Binary correction regularity and
axis/puncture basis compatibility remain to be derived before solving.

## Minimal acceptance table (declared before physical runs)

| Distinct obligation | Check and acceptance | Status |
|---|---|---|
| Seed conventions | Schwarzschild closed form; mass scaling; rotated/translated Kerr. Scaled field error<1e-11 at moderate parameters | pending |
| Derivatives and vacuum | Independent metric/K finite differences at3 steps; expected convergence before rounding; mass-normalized exterior H/M RMS<1e-7 | pending |
| Boost and causal slice | Generic spin/boost andGamma10; positive metric and slice margin at every evaluated point; E/P/J extrapolation error<1e-4 with two angular levels | pending |
| Linearization | One curved manufactured operator case and centered JVP difference sequence; relative error<1e-7 | pending |
| Moderate binary | Three increasing resolutions; independent fixed exterior near/bulk H/M RMS<1e-6 and max<1e-4 with decreasing errors; positivepsi; charge changes<1e-3 | pending |
| Consumer/horizons | Bound checkpoint and field/gradient roundtrip<1e-12; initial-time expansion RMS<1e-7 at3 angular orders, mass/spin changes<1e-4; modified-region enclosure | pending |
| QI compatibility | One unchanged saved generic QI fixture, default-path bitwise comparison on identical build/platform | pending |
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

## First implementation evidence

`src/HiSpID_trumpet.hpp` implements a standalone templated seed kernel using
the shared four-variable jets and geometry data type. It is not yet connected
to the public configuration, binary caches or checkpoint interface.
`tests/reference/kerr_trumpet_81c7d803.hpp` is the pinned source provider;
its BSD license is retained beside it asATHENAK_LICENSE.

The54-point algebra smoke check (spin0/.9/.99; speed0/.5/sqrt(.99)) has no
rejected slices, minimum sampled slice margin.0132698, maximum unboosted
reference scaled field difference2.203e-13 and native Hamiltonian2.214e-12.
This is limited sampled evidence, not a proof of global spacelike slicing.
The probe evaluates physical metric time derivatives after the boost.

The separate raw-tensor Cartesian observer gives finest-step H/M RMS:
static5.330e-10/3.889e-12, generic spin/boost2.107e-10/2.229e-12,
Gamma10 3.207e-8/9.905e-11, at three exterior points per case.
All meet the absolute1e-7 bound. Static/generic momentum and derivative
differences display the expected decrease. Hamiltonian differences already
encounter rounding; Gamma10 derivative differences do not yet show a clean
convergence interval. These results do NOT complete the derivative or full
seed acceptance rows. Retained files are under`validation/trumpet/`.
Next: resolve derivative conditioning with a suitable independent observer,
complete geometric/charge controls, then wire an ABI-safe seed-family option
through configuration, caches, sampling and versioned checkpoints.

## Integrated family dispatch and current evidence

The explicit family API now selects QI=0 or R0=M trumpet=1 for both active
holes. `HiSpID_Config` retains its original ABI and existing entry points
retain QI behavior. The family is carried by the context, host/device cache
builder, sampler and operator hook. Python accepts `seed_family='trumpet_r0_m'`
on creation or as configuration metadata; it rejects unsupported families
and conflicts with checkpoint metadata. Sampling-only contexts remain CPU.
GPU dispatch is wired but not compiled or qualified in this stage.

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
