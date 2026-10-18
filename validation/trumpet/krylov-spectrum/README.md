# One retained-state Arnoldi diagnostic

Both overintegration and two damped modal corrections fail their bounded
fine-grid trials. Before another preconditioner or Krylov-method change,
record the unrotated Arnoldi Hessenberg matrix from one original-system,
original-preconditioner replay of the retained240 state. The corrected memory
accounting permits128 vectors in the restart window within the64GiB budget;
this single bounded cycle supplies more spectral information than the old64
window. It is not a physical-validation study or a backend sweep.

PUNCTURE_KRYLOV_ARNOLDI is an opt-in path for a JSON-lines diagnostic. Before
Givens rotations, copy each small host Hessenberg column; append one matrix
per completed restart cycle. No native vector kernel, dot product, arithmetic
order, stopping tolerance, or default allocation changes when the variable is
unset. Existing manufactured nonsymmetric controls pass227 checks; diagnostic
off/on logs are byte-identical. Independent NumPy least-squares residuals of
all121 recorded small cycles agree with native recurrences to1.12e-17 scaled.
Default Bowen--York BiCGStab is untouched. Arnoldi matrices from augmented
LGMRES columns must not be interpreted as ordinary GMRES Ritz projections;
this retained-state diagnostic uses ordinary GMRES only.

The replay caps one Newton attempt at128Krylov, original relative target0.1,
original nonlinear tolerance1e-12. It records a changed diagnostic checkpoint
if a Newton update succeeds, avoiding repetition of a useful step. It does
not run physical observers or inherit physical acceptance. The recorder's
small matrices must pass an independent least-squares consistency check before
their spectra guide another change. Singular values/Ritz values of a finite
projection are diagnostics, not a proof of the full operator's spectrum.

Job59955993 completed the original-M128-step projection. No iterate changed;
true residual ratio0.54673415 reaches the previous long-run plateau within
one cycle. The projected norm is4.29e12 and its singular-value spread is about
5.3e15, with eigenvector condition about4.2e12. Tail singular values and Ritz
values are therefore sensitive to rounding; no precise eigenvalue-based
deflation count is inferred. An untruncated independent QR reconstruction
agrees with the native recurrence within7.9e-10 absolute relative residual.
Default SVD truncation discards directions in this ill-conditioned matrix and
must not be mistaken for a recorder error. The raw matrix, analysis driver
and diagnostic plot are retained. This motivates a left-preconditioning trial
with explicit confirmation of the unchanged original-system residual.
