# HiSpID backend: build, conventions and interface

This independent backend constructs non-conformally-flat binary Cauchy data
and solves four coupled CTT constraint equations. Its acceptance status and
supported cases are in [HISPID_STATUS.md](HISPID_STATUS.md); numerical evidence
is in `validation/results.json`. A successful native Newton status is not a
certificate of physical constraint accuracy.

## Isolated build and tests

Requirements: C99 and C++17 compilers, GSL with `gsl-config` on PATH, and
Python 3 with NumPy for the adapter and independent verifier. No SciPy,
Einstein Toolkit, AthenaK or Python-project installation is required.

```sh
cd /Users/hz0693/research/lazarus/.hispid-worktrees/TwoPuncturesC
make -j1 hispid
make -j1 test
make -j1 test-hispid PYTHON=python3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONPATH="$PWD/python:$PWD/validation:$PWD/examples"
python3 validation/run_validation.py --library "$PWD/build-hispid/libHiSpID.so" --stage seeds
python3 validation/run_validation.py --library "$PWD/build-hispid/libHiSpID.so" --stage moderate --far-radius 0 --label moderate_far0
python3 validation/check_covariance.py --library "$PWD/build-hispid/libHiSpID.so" --case moderate_far0
```

On the development host the explicitly used NumPy runtime is
`/Users/hz0693/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3`.
Use it as `PYTHON` if system Python lacks NumPy. `make clean-hispid` removes
only the new backend's build products. Existing BY targets and C sources are
unchanged; the new library reuses their mapping/spectral routines through
the static library, without changing the BY global parameters.

High-regime validation fails closed unless the saved seed, moderate and solved
coordinate covariance gates pass for the same library. The moderate command
above reproduces a current **failed** diagnostic; no regular-basis solved
binary is yet accepted. Historical preliminary passes used an incompatible
axis-nonregular representation. `--levels 24:12,40:20,56:28` means `(N_A,N_B,N_phi)` of
`(24,24,12)`, `(40,40,20)`, `(56,56,28)`. `--angular-ratio 2` doubles
only the polar size, records the actual nonsquare shape, and prolongs in
the identical mapped coordinates. Give separate case labels to different
refinement experiments. Results are saved incrementally;
failed iterates and physical samples remain in ignored `validation/raw/`.
Only one numerical job should run at a time on the shared host.

## Free data and signs

Units are G=c=1; signature is (-,+,+,+). Extrinsic curvature is
`Kij=(Di beta_j+Dj beta_i-dt gammaij)/(2 alpha)`.
Each hole input supplies rest Kerr mass m, rest Cartesian spin S,
lab velocity v and center. Require `|S|<m²` and `|v|<1`; mass zero disables
a hole. These are seed parameters. They are not binary horizon masses or
spins, and input m times Gamma is not a measured binary ADM energy.

The rest quasi-isotropic four-metric is pulled back with
`t0=Gamma(t-v·X)` and
`X0=X+(Gamma-1)v(v·X)/v²-Gamma v t`. Four-variable second-order automatic
jets include time derivatives after boosting. The normal uses a signed
lapse through the Einstein–Rosen throat. Rest extrinsic curvature uses the
axial Killing identity for beta=omega*l and analytically factors
`dR/dr=(r-c/r)/r`, `Delta=(r-c/r)²` from `(partial omega)/alpha`.
The boosted graph slice `t0=-v·X0` then gives its second fundamental form
without dividing by alpha. This is geometrically the same seed as the
four-metric ADM decomposition and remains finite at the QI throat; the
throat's physical constraints and K derivatives are tested. Exact punctures
remain excluded. `long double` follows the platform ABI (on this Apple ARM
host it has double precision), rather than providing arbitrary precision.
Near-extremal or extreme-boost binaries require additional validation.

For isolated seeds, `E=Gamma m`, `P=Gamma m v`, and
`J=Gamma S-Gamma² v(v·S)/(1+Gamma)`. The thesis's boosted-Kerr interface
uses lab spin; convert it by
`Srest=Slab/Gamma+Gamma v(v·Slab)/(1+Gamma)` before passing it here.

`conformal_choice=0` uses `(Sigma/r0²)^(1/4)`; choice 1 uses
`det(gamma_seed)^(1/12)`. These give the same isolated physical seed but
different binary free data. The final conformal metric generally has
nonunit determinant. With companion filter f and optional far filter F,

```
h = delta + sum(f F (h_seed-delta))
Psi = 1 + sum(F (psi_seed-1))
K = sum(f F K_seed)
Mraw = sum(A_seed),  A_seed = psi_seed² (Kij_seed-gammaij_seed K_seed/3)
M = Mraw - h trace_h(Mraw)/3
psi = Psi + u
Atilde = M + L_h b
gammaij = psi⁴ hij
Kij = psi^-2 Atildeij + gammaij K/3
```

Projection uses the superposed metric before raising indices; all derivatives
of the projection are included. Covariant A is never multiplied by f or F.
The source evaluation uses each exact seed's vacuum momentum identity,
`Div_seed A_seed=(2/3) psi_seed^6 grad_seed K_seed`. Metric and raised-tensor
differences are formed before differentiation, and the corresponding
connection difference is evaluated covariantly. The trace is computed from
the inverse-metric difference, using the exact seed trace-free identity.
This is algebraically the same superposed projection/divergence and avoids
cancelling singular isolated terms in floating-point arithmetic. Unboosted
QI Kerr has exactly K=0, which is imposed before forming A. Native direct
divergence comparisons and independent raw physical momentum checks verify
the reformulation.
Companion `f_h=1-exp[-(r_other/omega_h)^p]`; far
`F_h=exp[-(r_same/far_radius)^4]`. Nonpositive widths/radii disable them.
The equation attenuation g is the product of the published C-infinity
tanh(tan) windows on `[inner_min,inner_max]`; nonpositive maximum disables
that hole's window. Its endpoint evaluation saturates below double precision.

The solved equations are paper v3 Eq. (26):

```
FH = Delta_op u + g [Delta_h Psi - psi R_h/8 - psi^5 K²/12 + A²/(8 psi^7)]
FM = DeltaL_op b + g [Div_h M - (2/3) psi^6 grad_h K]
```

With `inner_flatten=1`, Eqs. (27)–(28) prescribe
`h_op=delta+g(h-delta)` and `Gamma_op=g Gamma_h`. This connection is not
the Levi–Civita connection of h_op. The implemented interior extension
forms covariant L with h_op/Gamma_op, raises with h_op inverse, then takes
the divergence with Gamma_op and its derivatives. Physical reconstruction
always uses actual h. At g=1 the operators are the standard curved CTT
operators. Modified-region residuals must be reported separately.

## Numerical representation

The binary map is centered at the midpoint and aligned with its separation
by a proper orthogonal frame. All sampled tensors and correction vectors
are returned in the input lab frame. The current experimental representation
stores real orthonormal Fourier modes of a regular auxiliary field P, with
Chebyshev interpolation in two mapped coordinates. With X the prolate radial
coordinate and R its polar angle, define `a=tanh(X/2)`, `t=a²`, `eta=cos(R)`
and `q=a sin(R)`. Each mode reconstructs as
`u_m-W_m=-2(1-t)q^r P_m`, and likewise for each vector component without W.
Here r=m through m=4; higher odd/even modes use r=3/4. This parity-preserving
cap guarantees Cartesian C2 limits for arbitrary bounded P without division
by tiny q^m. Smooth higher modes are represented by additional powers of q
in P. It does not establish spectral convergence for puncture free data.

The raw Chebyshev Gauss coordinates are `z,zeta` in (-1,1). With
`sigma=(1+z)/2`, the fixed maps are `t=.2 sigma/(1-.8 sigma)` and
`eta=tanh(2 zeta)/tanh(2)`. Both maps are analytic and monotone, with nonzero
endpoint Jacobians. Meridional derivatives use off-diagonal row differences
`sum_j D_ij(P_j-P_i)` and a zero-sum second-derivative diagonal, so constants
differentiate to exactly zero. This changes rounding, not the polynomial
interpolant. Derivatives include their full first/second chain rules
and analytic mode factors. Fourier second derivatives retain the cosine
Nyquist mode independently of the first derivative. Mode-indexed arrays are
never subjected to a second Fourier transform during coefficient conversion
or prolongation. The exact native identifier is
`modal_P_C2prolate_mapped_v2`; old nodal-V arrays are incompatible.
For resolution experiments, build into a fresh directory with, for example,
`make -j1 hispid HISPID_DIR="$PWD/build-hispid-focus05_k3" HISPID_MAP_FLAGS='-DHISPID_RADIAL_STRETCH=.05 -DHISPID_ANGULAR_STRETCH=3.'`.
The general maps are `t=lambda sigma/(1-(1-lambda)sigma)` and
`eta=tanh(kappa zeta)/tanh(kappa)`, with positive endpoint Jacobians.
Nondefault maps receive the identifier
`modal_P_C2prolate_map_v3_r<lambda>_k<kappa>`, encoding both double constants
at17 significant digits. `HiSpID_collocation_maps` also returns the constants.
Map changes need fresh physical convergence evidence; they cannot reuse
old modal arrays as an initial guess. Use a fresh build directory because
Make does not track changes to command-line compiler flags.
The infinity factor allows a vector 1/r term when correction momentum is
nonzero. Leading angular coefficients are determined by the elliptic solve;
finite-grid asymptotic charge convergence still requires independent checks.

An exact change of scalar variable removes the known far-filter shell:
`u=W-2(1-t)sum_m q^r P_m`,
`W=sum((1-F)(psi_seed-1))`. W and its first/second derivatives are computed by
automatic jets; W is O(r³) at a puncture and O(1/r) at infinity. Its puncture
extension is C2, with vanishing first and second derivatives, but need not be
analytic or C-infinity. This variable change alone does not establish
exponential spectral convergence. It is added to base
fields and omitted from JVP directions. The full `Delta_op W` is retained
even when g<1. This changes neither Eq. (26), free data, nor boundary
conditions. Finer validation grids can start from tensor-product
interpolation of a converged coarser auxiliary field; each level still solves
its own coupled equations to the declared tolerance. Old checkpoints made
before this change stored u rather than
u-W and require conversion. Checkpoints are tied to configuration and
library SHA; arbitrary unknown arrays are not portable between grids.

Newton uses an analytic Jacobian-vector product, restarted right-
preconditioned GMRES with two-pass orthogonalization and a damped line search.
A Fourier-mode FD preconditioner approximates the flat Laplacian of
the analytic mode factor, including its potential and mapped-coordinate chain
rules. It uses the azimuthal mean of `trace(h_inverse)/3`, the scalar Newton
potential and a 4/3 Laplacian approximation for vector rows. Its five-point
matrix is inverted by block-tridiagonal elimination with pivoted dense radial
blocks. Cosine/sine partners share factors; scalar and vector factors are
separate, and the three approximate vector rows share a factor. This removes
the slow long-wavelength convergence of the earlier ILU(0) approximation.
The full curved coupling remains in every residual/JVP. On the tested N80
moderate configuration, exact modal block elimination reduces the serial
solve from296.9s/896 Krylov iterations to31.3s/64, while retaining the same
failed physical residuals. More Krylov storage can help restart stagnation:
a constant-metric anisotropy control fails with restart32/max300 but its
two difficult m0 vector rows pass with restart128/max1024 at unchanged
1e-11 relative tolerance. Other failed rows remain unresolved. These
controls are available through `build-hispid/test_solver.x
--anisotropic-vector-only` and `--anisotropic-difficult-only`; they are
separate diagnostics from the default passing native suite. Compensated Fourier projection removes
the row mean from nonzero modes; exact Fourier rows have zero sum, so this
changes neither equations nor retained modes. Residual rows are multiplied
by `(sin(alpha) sin(beta))^6` as in the inherited solver. Both weighted and
raw conformal extrema are exposed. Neither replaces physical validation.

The required horizon-enclosure check concerns `g<1` and the optional modified
correction operators. The exponential companion/far profiles f/F modify
free conformal data and have noncompact tails; those tails cannot all fit
inside a finite horizon. At g=1 the full coupled equations still impose
vacuum wherever f/F differ from1. Record all profiles and independently
check the exterior physical constraints. Isolated Kerr throat radii only
screen candidate g windows; the solved binary finder must establish their
actual enclosure.
Allocation is rejected if a conservative cache/stencil/Krylov estimate exceeds
the per-context `memory_limit_mib` (default2048, accepted range16--8192),
including the additional dense block factors.
High-spin refinement beyond the initial grids explicitly uses4096 MiB,
or6144 MiB for160²×24;
only one context/job is active during those large solves. Match the native
header, adapter and library when using this experimental ABI. Sampling
differentiates the coefficient basis directly, including the off-grid cosine
Nyquist derivative. The inverse prolate map uses stable Cartesian distance
expressions near every axis segment. Exact-axis values and first gradients
use analytic m=0/m=1 limits, and transverse radii below1e-10 b use their
first-order Taylor extension. The foci themselves remain excluded. These
limits are tested against an independent Cartesian distance expression;
binary vacuum accuracy is a separate convergence gate and is still pending.

`HiSpID_create_sampler` loads the same saved unknowns without allocating the
collocation background, derivative workspace or Newton/Krylov work arrays.
Its allocation screen uses128 bytes per point, including coefficient-transform
scratch headroom. It rejects solve, residual, JVP and equation-sample calls.
`HiSpID_sample_with_derivatives` adds physical lab-frame metric first gradients,
stored as `dgamma[27*p+9*d+3*i+j]`; this includes the fixed analytic far
correction W and three frame rotations (derivative and both tensor indices).
Python exposes `Backend.create_sampler` and `Solution.sample_with_derivatives`.
The adapter checks the loaded symbol image and records its first-load SHA.
Compare distinct builds in separate processes: macOS can reuse an archived
image with the same embedded install name. The migration verifier does so
explicitly, and normal checkpoint replay still requires a matching build.
Archived libraries remain loadable for explicit migration comparisons; methods
requiring new symbols fail clearly if those symbols are absent.

The portable AthenaK interchange file is explicit-field text version1, written
by `examples/export_athenak.py`, with seventeen digit doubles, source-library
SHA, an acceptance label and a strict unknown count. It never serializes native
structure padding. The source SHA records evidence provenance; cross-platform
use additionally requires source/version agreement and validation of the new
build. The isolated AthenaK pgen imports physical gamma/K directly, checks the
ADM/Z4c round trip, and can call its finder at time zero. Direct native-geometry
finder checks and mesh-interpolation checks are reported separately.
Both map centers must be distinct, including a mass-zero inactive focus,
so every active puncture remains at a focus. The separate seed API does not
have this map restriction.

## C interface and Python use

`include/HiSpID.h` is a C ABI with owned opaque contexts. Initialize a config
with `HiSpID_default_config`, customize it, create, solve, sample and destroy.
`solve` returns 0 for weighted convergence, 1 for numerical failure and -1
for invalid context. Failed iterates remain sampleable. Other evaluations
return zero on success and a negative code on failure; consult
`HiSpID_last_error`. The context is not safe for concurrent calls because
it owns reusable spectral buffers. Independent contexts do not mutate BY
globals. Arrays of physical gamma, Kij, h and covariant Atilde have full
row-major 3×3 storage; XYZ is an array of triples. Point correction is
`[total u,bx,by,bz]`; mean curvature is the prescribed physical trace.

`get/set_unknowns` uses four interleaved values per collocation point:
`V0,Vx,Vy,Vz`, at index `4*(i+NA*(j+NB*k))+field` in the local map frame.
These hooks and the residual/JVP/operator hooks support reproducible
diagnostics, not a second physical-data interface. The adapter requires
an absolute library path:

```python
from hispid import Backend
from configs import moderate

backend = Backend('/absolute/worktree/build-hispid/libHiSpID.so')
config = moderate(backend, 40, 20)
config.far_radius = 0
with backend.create(config) as solution:
    diagnostic = solution.solve()
    physical = solution.sample([[1.1, .2, .3], [8., 1., -.7]])
    finite_radius_EPJ = solution.charges(200)
```

See `examples/integration.py` for a command-line example that saves arrays
and diagnostics without importing or modifying the Lazarus Python project.
ADM charges are general physical metric/K surface integrals on finite lab
spheres, returning `[E,Px,Py,Pz,Jx,Jy,Jz]`. Extrapolate across increasing
radii and check both polar and azimuthal quadrature. The integration polar
axis follows the prolate separation axis; normals, tensors and returned
charges remain in the lab frame. The energy flux uses the tested analytic
physical metric gradient, independently checked against Cartesian finite
differences. This reduces the sampler calls per node from thirteen to one.
J is about the supplied sphere center;
changing its origin from 0 to c gives `J_about_c=J_about_0-c cross P`.

`HiSpID_residual_scaling()` reports the actual row norm. The default remains
`sin6_alpha_beta`. The optional compile flag
`HISPID_EXPERIMENT_FLAGS=-DHISPID_INFINITY_EQUILIBRATION=1` multiplies it by
`(1-t)^-6`; this is an unsupported conditioning diagnostic. Its declared
exact-seed far-source floor fails at up to1.89e-13 versus1e-14. No binary
has been solved or accepted in that norm. Experimental flags require an
isolated fresh build directory, never a shared/default library overwrite.

## Independent validation and reference departures

`validation/physical.py` differentiates only sampled physical gamma/K with
fourth-order Cartesian differences and tensor-product mixed stencils.
It computes `H=R+K²-Kij K^ij`, all contravariant momentum components and
their physical norm. It reports g<1 and g=1 separately, near and bulk bins,
and normalized ratios; denominator floors are `1e-8` in the chosen mass
units. Adaptive high-regime steps are `min(.002,.001 distance_to_puncture)`
and are refined to assess the verifier floor.
These are RMS values over fixed off-grid point sets, not volume-integrated
L2 norms. The cancellation-normalized momentum ratio is uninformative for
maximal data: with K=0 its denominator is the already-cancelled divergence.
Absolute physical momentum norms define acceptance.

Equation/page derivations, primary sources, original implementation search,
source typos and independent review evidence are recorded in
[hispid-review.md](hispid-review.md). This backend targets the later v3
smooth source attenuation plus full trace projection. It does not silently
substitute the thesis's 2015 discontinuous tensor stuffing. The HS99UU input
uses thesis companion width .2 and QI scalar, but retains modern trace
projection; any charge disagreement must be reported. High-boost tables
give measured horizon masses and incomplete bare input parameters; local
benchmarks at the same boost are labeled accordingly.

There is no binary apparent-horizon finder. Seed coordinate horizon-radius
screening does not verify enclosure in the solved geometry. No horizon
mass/spin or horizon-contained exterior-vacuum claim is made.
