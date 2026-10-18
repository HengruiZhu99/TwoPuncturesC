# Left-preconditioning candidate rejected

The retained right-preconditioned128-step Arnoldi projection reaches the old
stagnation floor within one cycle. Its large upper couplings produce a singular
value spread near the double-precision limit (about5e15). Exact individual
Ritz values are too sensitive to infer a deflation dimension. See the raw
matrix, independent QR check and plot in ../krylov-spectrum/.

This opt-in HISPID_LEFT_GMRES changes the minimization metric instead of
changing the PDE or enlarging another restart window. For the original frozen
Jacobian A and modal inverse M, each inner cycle solves

    M^-1 A delta = M^-1 (rhs - A x)

with ordinary GMRES and an identity inner preconditioner. Update x += delta,
then explicitly recompute the residual in the original system. Success is
allowed only when its original Euclidean L2 norm is <= eta ||rhs||2. The
nonlinear sin^3/tau residual and1e-12 stopping bound remain unchanged, as do
all physical observers and horizon criteria. The inner iteration budget is
shared across cycles, not reset without limit. The same backend-independent
outer policy runs with native CPU or Kokkos vector primitives. Three/four
scratch vectors are included conservatively in the existing64GiB budget.
Default/OFF, QI and Bowen--York paths remain unchanged.

The local coupled control converges5Newton/75Krylov in0.587s. Frozen residual
and JVP are bitwise identical to the original operator; converged sampled
fields differ1.15e-13 (limit1e-10). This is an equation/solution control, not a
fine-grid performance claim. Native counters count actual Jacobian/modal
inverse applications, including original-residual confirmations. For this
side, execution controller operator_seconds includes the composed M^-1 A
operation, precondition_seconds refers to its identity action, and
linear_seconds sums inner cycles. Whole solve diagnostics.seconds is the
primary comparable solve time; do not interpret these inner phase timers
as exclusive raw A/M timings.

One CUDA small control gates one240-grid pilot from the retained failed state,
restart64,512total Krylov per Newton, eta=.1. Arnoldi matrices are retained
without another solve. Failed solves retain their checkpoint/counters but skip
redundant physical observers. If useful, the same frozen image/formulation is
assessed at224/240/256 before the finest AthenaK horizon check. No result from
the rejected Galerkin or damped-modal candidates enters that assessment.

The bounded production pilot59956254 completed:512 iterations, original
relative linear residual119652.944 (target0.1), solve114.070s. The original
residual confirmation correctly rejected this step. No Newton update was
applied; physical sampling was skipped and no finer study/AthenaK run follows.
The opt-in flag remains OFF. Exact result, Arnoldi cycles and log are preserved
in evidence240/. This failure redirects the investigation to regularity and
seed-family controls, rather than another Krylov parameter sweep.
