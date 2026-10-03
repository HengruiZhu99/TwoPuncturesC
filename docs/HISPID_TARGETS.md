# Revised validation and AthenaK targets

The October 2026 follow-up requests a Kokkos port and comprehensive matched
performance report before separate aligned seed chi=.99 and inward head-on
input Gamma=10 studies. This supersedes the earlier GPU deferral below.
The new quantitative plan is `validation/extreme_kokkos_plan.json`; every
physical result is currently pending. Earlier .95/.885 evidence remains
bound to those cases and cannot qualify the new targets.

For the boost study, start with equal rest masses .5 and d=25, then use the
measured component areas to report and, if needed, calibrate d/Mirr about50.
At Gamma10 the isolated seed energies sum to10, so d/Eseed is only2.5.
Measure the global ADM energy and both component masses. This is separate
from thesis Section4.5/Table4.3 (printed pp.117--119, PDF pp.142--144), whose
d/MADM=100--400 cases were chosen to approximate isolated holes. Retain all
changed-d cases, bounded component/common finder attempts and failed gates;
a failed common search does not prove absence.

Portable GPU-produced data will be checked against a pure CPU sampler in
separate processes, with both library hashes, exact basis/maps, checkpoint
hash and physical field/derivative/horizon-point witnesses. AthenaK's budget
reader, domain and finder-role changes are pending validation. The producer
budget is preserved; a consumer migration or ADM/Z4c round trip cannot
inherit physical acceptance. No evolution steps are authorized here.

The user's revised targets are seed rest spin |S|/m²=0.95 and lab coordinate
speed |v|=0.885 (Gamma=2.147807737), followed by initial-data import and horizon
checks in AthenaK. These inputs are not measured horizon spins or momenta.
Test separate spinning and boosted binaries and at least one combined case;
separate tests do not establish their simultaneous range. Generic moderate
unequal-mass cases and coordinate covariance remain required controls.

Use the same physical residual limits as the original strong gate: near/bulk
RMS1e-6, maxima1e-4, positive metric/conformal factor and weighted nonlinear
residual1e-14. Preserve the original all-norms monotonic flags separately.
For the new target runs, resolution improvement is required whenever a norm
is resolved above the independently calibrated verifier floor. A plateau
below that floor is reported as unresolved, not as a measured zero. Establish
the floor using exact vacuum controls on the identical point sets, multiple
Cartesian stencil sizes and recorded grid differences before making a
noise-qualified convergence claim. Charge angular quadrature and radial
extrapolation need separate refinement; use absolute uncertainty1e-5 as the
target, and retain the original preliminary0.5% stability flag.

The combined Kerr seed horizon screen is
r_h=m sqrt(1-chi²)/2 and minimum lab radial distance r_h/Gamma. For m=.5,
chi=.95,v=.885 this minimum is.0363452. Inner g radii must be scaled from
this Kerr value, not the spinless Schwarzschild radius. Window changes are
distinct configurations and never reinterpret an earlier failed benchmark.
This screening does not certify binary horizon enclosure.

AthenaK integration stays in
`/Users/hz0693/research/lazarus/.hispid-worktrees/AthenaK`, branch
`codex/hispid-pgen`, based on PR790 head
`22baa243970fa1880b2bbc48e88a590069d55e47` of
`HengruiZhu99/athenak:project/z4c_overhaul`. Original checkouts, shared builds,
and simulations remain untouched. Initial-time finder calls must not require
an evolution step. Validate physical gamma/K import and ADM/Z4c round trips,
then exact Schwarzschild/Kerr/boosted horizon controls and the solved binaries.

Finder success requires expansion-residual and angular-refinement checks,
not only irreducible-mass stabilization. AthenaK's current `hrms` column is
mean square expansion; report its square root explicitly. Require complete
angular coverage. A sufficient attenuation-enclosure condition is
|c_AH-c_hole|+r_gmax < min R_AH with a numerical margin and angular refinement.
Any mesh puncture floor/fill needs its own enclosure bound. The current
coordinate-rotation spin estimate is not an approximate-Killing-vector spin
for generic boosted/distorted surfaces.

Prioritize single-thread correctness. Profile elapsed time in geometry,
spectral differentiation, JVP, preconditioner and Krylov reductions before
choosing optimizations. Perlmutter access is optional, through the user's
existing SSH control master; numerical jobs use allocated compute nodes,
not login nodes. OpenMP/MPI may be added when measured runtime is a bottleneck.
GPU work is deferred unless CPU approaches are insufficient. All new jobs,
builds, data and logs stay in dedicated task directories.
