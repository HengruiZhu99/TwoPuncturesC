# Regularity and seed-family control

The user explicitly permits another conformal-metric family if trumpet
regularity contributes to the moderate-spin problem. This is one bounded
comparison using the existing quasi-isotropic (QI) Kerr implementation,
not another backend/performance matrix.

For the unmodified, nonspinning R0=M trumpet, linearizing the scalar CTT
equation with fixed free data gives the leading radial operator

    Laplacian(u) - 11 u/(4 r^2).

Its bounded homogeneous spherical-harmonic branch has exponent

    p_l = -1/2 + sqrt(l(l+1)+3),   p_0 = 1.2320508075688772.

This permits fractional powers even after the leading seed singularity is
factored out. It does not establish the coupled spinning/boosted indicial
spectrum or prove the observed binary error comes from this exponent.
In the implemented g=0, inner_flatten=1 plateau, the correction equations
instead reduce to the flat Laplace/Lame system. Thus this fractional local
potential is absent there. The C-infinity compact transition is nonanalytic,
and can still limit global spectral convergence. Its sharpness, the capped
axis representation, and ill-conditioned modal inverse remain separate issues.

The frozen control keeps the original generic m=0.6/0.4 binary, spins,
velocities, separation6, exact seed scalar choice, no far attenuation, same
maps/basis and original modal GMRES. Only the seed family and corresponding
isolated horizon-radius scale of the interior window change: fractions
0.15/0.5 of the QI throat radius, which is half the trumpet coordinate radius.
The QI solution is physically distinct; it is not expected to reproduce the
trumpet correction, charges or horizon masses exactly. No trumpet warm start
or acceptance is imported.

Use the retained original-M CUDA image (hash in plan.json), host geometry,
and CPU sampler/AthenaK consumer. Predeclared grids are160/192/224 with twice
as many polar points and16 azimuthal points. The160 pilot must meet the
unchanged physical bounds before launching the other two. Near-hole sample
radii follow the selected family's horizon scale, so direct QI/trumpet norm
ratios are not a controlled slicing-only comparison. Independent bounds,
bulk points and within-family convergence are meaningful comparisons.

Even a successful QI study would validate a separate fallback formulation,
not repair or qualify the trumpet family. A failed QI control will also be
retained; no family is accepted from its internal weighted norm alone.
Multigrid is deferred until there is evidence that a coarse-error component
is the bottleneck, rather than the measured grid-scale error/conditioning.


The narrow QI160 control59956570 internally converged in8Newton/46Krylov,
42.386s, but bulk momentum RMS3.80871e-6 fails the1e-6 bound. Near H/M are
2.50060e-7/5.67512e-7, bulk H8.73073e-9. No finer narrow-QI run followed.
The exact result and physical arrays are in evidence160/.

A separate wider-QI candidate changes only the outer transition fraction to
0.85 (plateau0.15 unchanged). Its48near/18bulk observer positions are bitwise
identical to the narrow QI observer; no points are moved outwards. The frozen
wide-plan.json again specifies160/192/224 and unchanged physical/horizon gates.
The coarse160 bounds pass, with near H/M3.13058e-7/9.14461e-8 and bulk
H/M6.30273e-9/1.04547e-7. This36.4-fold momentum improvement is a transition
comparison within QI, not evidence of completed three-grid qualification.

The offline retained_decay.png and retained-decay.json use existing160 fields
only. The m=0 factored coefficients have slow meridional tails in both families.
Degree8 polynomial windows localize the trumpet wide-case radial vector tail
mainly to compactified infinity: suppressing that endpoint reduces the selected
80--143 degree-band norm to0.00406 of its original value; suppressing the inner
radial axis leaves1.003. The narrow QI tail instead lies mainly near the inner
radial axis (which includes the foci and the segment between holes). This does
not identify exact logarithmic powers, prove coupled indicial exponents, or
compare the families with identical physical interiors. No filtered field is
used as solved data. Necessary low-m Cartesian corner-compatibility defects
are also recorded as representation diagnostics, without acceptance claims.


The192 grid meets absolute physical bounds, but bulk H/M rise slightly to
8.62285e-9/1.09038e-7. The original strict-decrease gate is therefore not
satisfied. The224 grid is retained to finish the already frozen study.
The user's explicit finest-data import/horizon request will still be carried
out diagnostically if all three grids meet absolute physical bounds. The
original failed trend flag stays failed; even successful horizons cannot
promote the binary to the original full qualification. This avoids repeatedly
postponing the independently requested AthenaK check. All finder, migration,
angular refinement and enclosure thresholds remain unchanged.


All three wider-QI grids pass absolute physical bounds, but bulk errors increase.
The original three-grid assessment therefore remains passed=false. The
supplemental combined66-point physical RMS decreases2.835e-7 ->1.406e-7
->1.033e-7; it does not replace the failed regional trend gate.

Finest224 migration59957080 passes at2454 points with bitwise-identical fields
and metric gradients. The first finder launch stops before any numerical
finder work because the Python proof checker omitted the current v5 identifier.
AthenaK commit e54102e5 admits explicit v4/v5 identifiers, keeps exact basis/map
and image bindings, and rejects unknown versions and invalid matching maps.
The v2 attempt59957101 reuses the passed migration and succeeds at orders
8/12/16 plus independent quadrature. Final MH are0.615673713/0.416086392,
Mirr0.604595771/0.409618329, coordinate spins0.370895278/0.345813441.
Final expansion RMS2.041e-8/5.626e-9; import error4.083e-16. Buffered continuous
enclosure margins0.03205/0.01867 are positive. All consecutive angular/quadrature mass/spin changes are
<2.01e-12/<2.79e-11. Full raw surfaces, logs, inputs, migration witnesses and
the failed first attempt are retained in horizons/. These measurements are
diagnostic, not a successful strict three-grid or trumpet qualification.
No multigrid or further backend benchmark is launched.
