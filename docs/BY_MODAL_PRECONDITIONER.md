# Scalar BY modal preconditioner

`TP_preconditioner=1` changes only the inverse used by BiCGStab. The original
spectral nonlinear residual, spectral Jacobian-vector product, update
`V -= dv`, and native stopping rules remain available unchanged. The default
is the inherited line-sweep preconditioner (`TP_preconditioner=0`). The new
inverse is serial and requires even nphi>=4, as documented in `TP_Modal.h`.

## Analytic meridional stencil

The original solver stores V and the physical correction is u=dV, with d=A−1.
At native Gauss nodes A=−cos(alpha), B=−cos(beta), define

```
a=(A+1)/2, X=2 atanh(a), R=pi/2+2 atan(B)
h=1-a², hX=-a*h, k=(1+B²)/2, kR=B*k
H=b² (sinh²X+sin²R), rho=b sinhX sinR
```

The physical Euclidean Laplacian in these orthogonal prolate coordinates is

```
Delta = (dXX+coth(X)dX+dRR+cot(R)dR)/H + dphiphi/rho².
```

The vacuum BY scalar Jacobian on u is `Delta - C`, with
`C=(7/8) BY_KK(x,y,z)/psi^8`. Thus its action on dV has coefficients

```
cAA = d*h²/H
cA  = [2h² + d*(hX+h*coth(X))]/H
cBB = d*k²/H
cB  = d*(kR+k*cot(R))/H
cphi= d/rho²
c0  = (hX+h*coth(X))/H - d*C.
```

For example the alpha coefficients are
`q2=cAA/sin²(alpha)` and
`q1=cA/sin(alpha)-cAA*cos(alpha)/sin³(alpha)`.
The beta coefficients follow the same chain rule. With row weight
`w=(sin(alpha)sin(beta))³` and steps ha=pi/na, hb=pi/nb,
neighbor coefficients are `w*(q2/ha² ± q1/(2ha))`, and similarly in beta.
The diagonal subtracts the two neighbors in each direction and adds w*c0.
Ghost neighbors follow the original `Index` reflection: each outward i/j
neighbor folds into the SAME row diagonal, with no phi shift or mode parity.
Mixed derivatives cancel analytically; the original Cartesian transformation
can leave roundoff mixed coefficients, so this is not a bitwise copy of JFD.

Only C depends on phi. Average the POINTWISE C at fixed i,j, using current
physical `u->d0`, current bare masses, and the unsmoothed `BY_KKofxyz` used by
`LinEquations`. Averaging psi or KK separately would change the approximation.
After Fourier transformation, the original second-order FD phi eigenvalue is

```
lambda_m = -4 sin²(pi*m/nphi)/(2pi/nphi)².
```

Add `w*cphi*lambda_m` to the modal diagonal. This is the FD preconditioner
symbol, not the full spectral −m² symbol. Cosine/sine partners share a real
meridional matrix; the constant and Nyquist modes have one partner each.
The real transform is orthonormal, with normalization1/sqrt(nphi) for those
two modes and sqrt(2/nphi) otherwise. Nonzero modes subtract a compensated
azimuthal mean to reduce cancellation; the inverse maps back to NODAL V.

## Block inverse and failure semantics

Each frequency has a block-tridiagonal meridional matrix (polar j is the block
index; each block has na radial entries). With blocks L_j,D_j,U_j, elimination
uses `S_j=D_j-L_j*T_(j-1)` and `T_j=S_j^-1 U_j`. GSL pivoted LU factors S_j.
Forward/back substitution solves the whole meridional system. All rows are
first equilibrated by their largest coefficient; the RHS receives the same
scaling. Dense lower/transfer blocks retain Schur fill. Factors are rebuilt
for every Newton linear solve, including target-mass adjustment; no stale
nonlinear potential is reused.

Three dense banks cost `3*(nphi/2+1)*nb*na²*sizeof(double)` plus transforms,
permutations and scratch space. Checked sizes bound total indices to INT_MAX
and private allocation to4GiB. Invalid/nonfinite inputs, singular pivots,
nonfinite factors or solves, and allocation failure reject the requested
modal solve without silently selecting line sweeps. Context scratch and the
inherited process-global BY parameters/GSL handler make this path serial and
non-reentrant. A nonzero pivot alone does not certify conditioning.

For opt-in modal or RHS-relative forcing, the linear gate recomputes TRUE
`||F-Jdv||₂` with the ORIGINAL spectral JVP. A failed gate rejects the Newton
step before changing V; target-mass adjustment stops before changing masses.
The retained context remains available for diagnosis, with diagnostics status1
and explicit linear/modal failure counters. Selectors must be0/1, and
`TP_linear_rtol` must be finite with0<rtol<1.

## Independent verification

`tests/test_by_modal.c` builds a FULL PHYSICAL-GRID reference by translating,
averaging and symmetrizing original `SetMatrix_JFD` rows around phi. It does
not call the analytic/modal assembly to build its reference. The legacy-JFD
constructor and analytic constructor must both invert this matrix on
manufactured meridional/Fourier fields, including every partner and Nyquist,
rectangular grids, reflected faces/corners and updated masses/corrections.
It also checks nonfinite/singular/malformed inputs, reuse after rejected RHS,
checked large dimensions and rejected target-mass steps. Current bounds are
2e-9 on solution error and2e-12 on normalized matrix residual; observed maxima
are3.342e-14 and1.772e-13. No spectral solve is replaced by this FD oracle.

The full40×80×16 benchmark and original standalone controls separately require
recomputed nonlinear residual<=1e-12, no opt-in failures, final fields and
charges within the predeclared1e-10 bounds, and independent Cartesian H/M
agreement within FD-floor bounds. Default full-grid states and small original
standalone cases remain bit-identical; modal core values, coefficients and physical samples are numerically equivalent.
Retained native near-puncture derivative arrays amplify roundoff/path differences
and are diagnostics outside the1e-10 core equivalence gate.
See `by_modal_acceptance.json`, `common_stopping_moderate_40.json`,
`by_modal_verification_v2.json` and `HISPID_PERFORMANCE.md` for provenance and
measurements. Computational stopping comparisons do not confer binary physical
validation or equal physical accuracy between BY and HiSpID.
