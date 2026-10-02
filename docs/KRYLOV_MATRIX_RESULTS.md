# Four-way binary Krylov comparison

All four methods reach the same computational stopping convention on the
moderate unequal-mass binary: 40×80×16, zero guess, one CPU thread, cubic row
weighting, outer maximum1e-12 and **true RHS-relative linear L2 <=1e-3**.
Each system uses its own fixed modal preconditioner, unchanged between methods.
GMRES restart64 and both caps2000. Two fresh-process interleaved repetitions
on Apple M5 Pro; solve seconds exclude construction and post-timing checks.

| System | Krylov | Newton | Iterations | Spectral JVPs | M applies | Median solve | Peak RSS |
|---|---|---:|---:|---:|---:|---:|---:|
| HiSpID | GMRES |4|42|50|42|3.521 s|304.35 MiB|
| HiSpID | BiCGStab |4|23|54|46|3.699 s|276.13 MiB|
| BY | GMRES |4|16|24|16|4.048 s|91.38 MiB|
| BY | BiCGStab |4|10|25|17|3.977 s|89.63 MiB|

BiCGStab iterations generally perform two JVP/M actions; comparing iteration
counts alone exaggerates its advantage. Hi BiCGStab saves **9.27% RSS** and
is **5.04% slower** in the final pair. BY GMRES uses **1.96% more RSS**;
there is no convincing time advantage for either BY algorithm. BY GMRES ranged
3.830–4.265 s; earlier retained experiment pairs gave a small advantage in the
opposite direction. These short trials are not a scaling curve.

The shared BiCGStab workspace also removes the duplicate derivative workspace
in inherited BY. Its native-default peak RSS is83,607,552 bytes versus the
archived comparison's91,668,480 bytes, **8.79% lower** in those fresh workers.
That memory comparison is one default run here versus the archived reference;
timing is sensitive to local load and no additional default speedup is claimed.

## Preservation and numerical limits

Both archived default solutions pass **complete bitwise state preservation**,
counts, and BY's full-precision iteration trace. Hi keeps GMRES/adaptive
forcing/sixth-power row weights; BY keeps BiCGStab/original forcing/line sweeps.
The initial extraction experiment failed this gate because inlined reductions
changed compiled rounding. The final BY adapter calls its original external
norm/dot routines; GMRES retains the original direct reductions.

BY GMRES versus BiCGStab passes the predeclared core values/coefficients,
ordinary physical metric/K/lapse/psi, charges and independent fourth-order
Cartesian constraint differences. Near-puncture native derivative arrays
remain finite diagnostics outside that core equivalence gate, as in the prior
modal study.

Hi BiCGStab passes native convergence, true linear residual histories,
sampled physical fields, charges and independent physical constraint checks.
It **fails** the predeclared1e-10 full auxiliary-P/coefficient equivalence bound:
raw-P scaled difference is about2.50e-5 and Chebyshev coefficient difference
about6.49e-8. The largest difference is cosine m6, vector component2, at the
first radial/angular node. The reconstruction multiplier there is about3.8e-17.
Across all original collocation nodes, reconstructed physical correction
values differ by at most about3.4e-13. This is evidence of weak observability
of the parameterization near axes/foci, not an exact-nullspace proof or a
certificate for global gradients/Hessians. The strict gate remains failed.
Hi BiCGStab is available as an explicit option; GMRES remains the default.

A small native sixth-power-weighted Hi fixture also failed raw-P1e-10 agreement
at common inner1e-3, and tighter inner1e-6 did not resolve it. Those failed logs
are retained. The public interface unit test checks sampled physical fields;
the separate binary report continues to enforce and report the failed saved
representation gate. No tolerance bound has been waived in that report.

## Reproducibility

- Criteria: `validation/krylov_matrix_acceptance.json`.
- Final matrix: `validation/krylov_matrix_moderate_40_final.json`.
- Earlier rejected extraction runs: `validation/krylov_matrix_moderate_40.json`
  and `validation/krylov_matrix_moderate_40_v5.json`.
- Values-only reconstruction diagnostic: `validation/krylov_representation_final.json`.
- Original BY fixtures: `validation/krylov_backend_fixtures.json`.
- Final suite/failure fingerprints: `validation/krylov_backend_verification.json`.
- Raw states, stdout and immutable timing images remain in their recorded
  `validation/raw/...` and `build-krylov...` paths in the isolated worktree.

The matrix driver returns failure because its strict Hi representation gate
fails, even though all four numerical solves converge and both defaults are
preserved. No physical high-spin/high-boost binary acceptance or horizon
qualification is inherited from these performance tests.
