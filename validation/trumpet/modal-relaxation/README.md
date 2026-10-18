# Grid-scale preconditioner correction candidate (not yet fine-grid qualified)

The retained failed240 residual has93.60% radial/78.56% polar energy in the
last eight Chebyshev modes; it is not a demonstrated coarse-grid error. The
existing preconditioner already solves its global modal FD approximation.
A second stationary residual correction therefore tests a high-frequency
smoother before adding a multigrid hierarchy.

The opt-in HISPID_DAMPED_MODAL_RELAXATION applies two damped corrections of
one frozen Jacobian system, z0=.4 M^-1 r and
z1=z0+.4 M^-1(r-A z0). It changes only the right preconditioner. The shared
backend-independent helper is instantiated with CPU and Kokkos operations;
all GPU workspace is allocated once per linear solve. The old path and QI
family retain one undamped application. Native work counters include both
FD solves and the extra Jacobian application. A native description query
records the resolved preconditioner in Python diagnostics. Trace mode emits
one recurrence-progress line every100 iterations; the controller still checks
the true residual before success.

This is a fixed linear polynomial preconditioner, not multigrid and not a
change to the nonlinear equations. The .4 damping is motivated by the
high-frequency spectral/FD Laplacian discrepancy, not fitted to a binary
physical observer. Curved coupled systems have no asserted spectral bound;
a retained fine-grid state must establish whether this is useful.

Small CPU control:5Newton/33outer Krylov,0.454solve seconds. Both frozen
residual/JVP are bitwise identical to the old operator, and converged fields
differ by1.68e-13 against the1e-10 control limit. QI residual/JVP remain
bitwise exact; centered trumpet JVP errors1.34e-11/8.12e-12. These controls
are not physical binary validation. Two sweeps cost two FD inverses and one
extra Jacobian per outer preconditioner call; iteration count alone is not a
speedup measure.

Before production use: one CUDA small-field comparison, then one bounded
240-grid pilot from the actual retained failed state, original formulation
and original observer. No backend/method sweep. If useful, validate a single
frozen image at224/240/256 and measure the finest horizon in AthenaK. The
Galerkin pilot is evaluated separately; its results cannot be combined into
this different discretization's study.

CUDA build59955507 stops at compilation: NVCC disallows an extended device
lambda nested in a generic host lambda. The scale operation now calls the
existing typed Kokkos vector helper. No numerical control or binary solve ran
in that job. Its build-v1 log is retained, and the retry writes build-v2.
The attempted overlapping CPU build on59955120 never launched because its
required ports were in use; no compute work ran in that step.

CUDA retry59955583 passes5Newton/33outer Krylov, 0.063434solve seconds; CPU/CUDA field scaled difference 3.44984e-13. Producer image7abe943bb14d35b49dc79599431abf4cd7666eb072ed317387009e49c3fb9c03 is frozen for the bounded240 pilot and conditional224/240/256 study.
