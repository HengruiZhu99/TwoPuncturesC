# Physical ADM sampling interface

This fork adds arbitrary-point sampling for a solved vacuum Bowen–York
puncture configuration. The underlying Hamiltonian solve and spectral
coordinates remain those of the upstream implementation at
`ec563aeb672235b9443c330f9cde65f7246e8ea4`.

```c
int TwoPunctures_sample_points(
    ini_data *data, int npoints, const double *xyz,
    double *lapse, double *psi_full, double *gamma6, double *K6);

int TwoPunctures_diagnostics(
    ini_data *data, double *residual_linf,
    double *adm_mass, double *puncture_masses2);
```

`xyz` has row-major shape `[npoints,3]`. `lapse` and `psi_full` have length
`npoints`. Both tensor outputs have row-major shape `[npoints,6]`, with
components **xx, xy, xz, yy, yz, zz**. Callers allocate all output buffers;
buffers must not overlap. Zero points is valid and needs no buffers. Return
codes are 0 for success, −1 for an invalid argument/context, and −2 for an
unsmoothed puncture or unsupported parameters. On failure, earlier output
rows may already have been written. `multiply_old_lapse=1` is unsupported
because this API has no input lapse field to multiply.

The full solved conformal factor is

$$
\psi=1+\frac{m_+}{2r_+}+\frac{m_-}{2r_-}+u.
$$

The physical covariant spatial tensors returned are
$\gamma_{ij}=\psi^4\delta_{ij}$ and
$K_{ij}=\psi^{-2}\bar A_{ij}$, with trace $K=0$ and the two-puncture
Bowen–York $\bar A_{ij}$. The sampling API returns physical tensors
independently of `conformal_state`.

The older Cartesian interface uses a different storage convention: for
`conformal_state>=1`, its `psi` is only
$p=1+m_+/(2r_+)+m_-/(2r_-)$, its `gij` is
$(\psi/p)^4\delta_{ij}$, and its conformal derivative arrays are derivatives
of the static factor normalized by $p$. Consequently physical
$\gamma_{ij}=p^4 g_{ij}$; those derivative arrays cannot be used as derivatives
of the full solved $\psi$. This convention is retained for existing users.

The lapse follows the legacy `initial_lapse` option. In particular `psin`
(enum value 2, the default) returns $p^n$, default $n=-2$, using the static
factor. `antisymmetric` returns $(2-p)/p$, `averaged` returns $1/p$, and
`brownsville` returns $2/(1+p^n)$. Lapse is a gauge choice and no shift is
specified by this API. A downstream ADM time derivative needs its own
explicit lapse and shift prescription.

## Parameters, units, and lifecycle

The punctures are at $(+b,0,0)$ and $(-b,0,0)$ before the existing
`center_offset*` translation and optional `swap_xz` interchange. Coordinate
separation is $2b$, not $b$. Input `par_m_plus/minus` are bare puncture
masses. With `give_bare_mass=0`, `target_M_plus/minus` instead request the
ADM masses at the two internal asymptotic ends. These are distinct from
apparent-horizon masses and from the ADM mass at the outer infinity.

`par_P_plus1..3`, `par_P_minus1..3` are Bowen–York momentum components;
`par_S_plus1..3`, `par_S_minus1..3` are angular-momentum components.
In geometric units, masses and momenta have length dimension, spins have
length squared, and $K_{ij}$ has inverse-length dimension. Momentum is not a
coordinate velocity, and spin input is not a dimensionless horizon spin.

Parameters are process-global and the library is not thread-safe. Use
exactly one live solve per process, serialized by the caller:

1. Set defaults and the physical and spectral parameters.
2. Call `TwoPunctures_make_initial_data` once.
3. Sample and read diagnostics without changing solved physical/spectral
   parameters. Affine frame and storage controls are interpolation options.
4. Call `TwoPunctures_finalise` on that context before setting up another.

Every solve now owns fresh allocations. A second live solve returns NULL;
`finalise(NULL)` has no effect. After finalization a new solve may use
different parameters and resolution. Sampling verifies that the current
spectral dimensions equal those used to allocate its context. A finalized
context must never be reused.

`TwoPunctures_diagnostics` returns the measured maximum absolute final
collocation equation residual, outer ADM mass, and internal-end masses in
plus/minus order. The residual is the actual stored equation residual, not
the requested Newton tolerance. This does not replace Cartesian exterior
Hamiltonian/momentum constraint convergence tests or horizon measurements.

With all momenta and spins zero, $\bar A_{ij}=0$ and $u=0$: these are exact
time-symmetric Brill–Lindquist data, with $K_{ij}=0$ and outer ADM mass
$m_++m_-$. A sufficiently close binary can have a common apparent horizon;
coordinate separation alone does not establish that. Constraint-satisfying
data also do not by themselves guarantee a valid linear Kerr perturbation.

## Native regression

Build with GSL available, then run `make test`. Shared and static libraries
exclude the standalone executable's `main` object, and objects depend on
the header. The static archive is recreated to remove any stale `main`
member from an older build.

The native regression performs three sequential solves: exact equal-mass
rest data, all-component nonzero boosts and spins, then different unequal
masses, closer separation, and smaller spectral dimensions. It recomputes
`F_of_v` into independent storage and checks the measured residual. It also
checks exact Brill–Lindquist fields and charges, an independently implemented
analytic Bowen–York tensor, nonzero $u$ in the full conformal factor, and
agreement with all four legacy conformal storage states and x/z interchange.

For the small generic-spin test, the recomputed residual is about
$9.0\times10^{-15}$. That is a solver regression at collocation points;
spectral truncation accuracy must be assessed separately for a scientific
configuration.

Primary formulation:
[Ansorg, Brügmann, and Tichy, gr-qc/0404056v2](https://arxiv.org/pdf/gr-qc/0404056v2).
