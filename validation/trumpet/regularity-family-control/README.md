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
