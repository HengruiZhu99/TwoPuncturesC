# Independent HiSpID formulation and validation review

Review started 2026-10-01. The reviewer did not edit implementation files or the original checkouts. This document is a review aid, not evidence that a solver has passed the gates below.

## Sources and availability

- [Ruchlin et al., arXiv:1410.8607v3](https://arxiv.org/abs/1410.8607v3), PDF pages 3–9 (formulation/numerical solve), 10–17 (initial-data cases/convergence), 20–22 (boosted evolution examples), 22–23 (ADM charges). Relevant sections were read from the locally downloaded full PDF extraction, with equations checked against the [arXiv HTML](https://arxiv.org/html/1410.8607v3).
- [Ruchlin dissertation, RIT 8797](https://repository.rit.edu/theses/8797/), supplied by the user as the complete 209-page PDF after the current and legacy repository download endpoints returned HTTP 403. Printed pp. 48–57, 64–70, 74–81, 97–106, and 128–140 were read in full; the relevant published-case sections on printed pp. 83–86, 107–119, and 141 were also checked. Printed-page/PDF-page offset is +25. The supplied source is stored in ignored `reference/thesis.pdf` and `reference/thesis.txt`. The prior access limitation is resolved.
- [Healy et al., arXiv:1506.06153v2](https://arxiv.org/html/1506.06153v2), sections II–III: supplementary high-boost construction and initial-data convergence. This is the primary paper cited as reference 54 in 1410.8607v3. Its earlier attenuation prescription differs from the later v3 formulation, so it must not silently replace Eq. (26).
- Original implementation search: GitHub repository queries for HiSpID, Ruchlin, and TwoPunctures; exact physics-term web searches; and public EinsteinInitialData branches. No publicly accessible original HiSpID backend was found. Unrelated botany repositories named hispid were excluded.
- The public [EinsteinInitialData repository](https://bitbucket.org/einsteintoolkit/einsteininitialdata/) has a `momentum_constraint` branch, commit `92c0f0a3067b8bb3a9826037768a632c16d6648e` (Roland Haas, 2010). Inspection of `TwoPunctures/src/Equations.c` and `TwoPunctures.c` shows conformally flat BY plus matter sources, with four Poisson potentials, not curved HiSpID. It is not an equivalent solver. The shallow reference copy lives in ignored `reference/et-momentum-constraint`; no code was copied from it.

## Conventions and geometric derivations

Signature is (-,+,+,+), geometric units G=c=1, and

\[
K_{ij}=-\tfrac12\mathcal L_n\gamma_{ij}
       ={D_i\beta_j+D_j\beta_i-\partial_t\gamma_{ij}\over2\alpha}.
\]

This is paper Eq. (8), PDF p. 4. It fixes the sign of the momentum charge. Given the four-metric and its first derivatives, an independent algebraic route is

\[
K_{ij}=-\alpha\,{}^{(4)}\Gamma^0{}_{ij},\qquad
\alpha=\pm(-g^{00})^{-1/2},\quad \beta_i=g_{0i}.
\]

The lapse sign chooses the normal orientation and must be carried through the Einstein–Rosen throat. For Schwarzschild, the smooth lapse is explicitly `(1-m/(2r))/(1+m/(2r))`, not its absolute value. For QI Kerr let `c=(m²-a²)/4`, `R=r+m+c/r`, `Sigma=R²+a² cos²theta`, `Delta=R²-2mR+a²`, and `A=(R²+a²)²-a² Delta sin²theta`. A signed lapse can be written without an absolute-value ambiguity as

\[
\alpha_0=(r-c/r)\sqrt{\Sigma/A}.
\]

Its zero occurs at `r_H=sqrt(c)`. Using positive lapse on both sides flips the inner-sheet K and introduces a discontinuity, even though the outer seed may pass vacuum checks. The exact throat has a coordinate degeneracy; a zero-lapse formula requiring division must be evaluated by its limit or avoided. A random grid missing the exact throat is not a verification of throat handling.

For a boost with input rest-frame mass m and seed momentum P, use `v=P/sqrt(m²+|P|²)`, `Gamma=sqrt(1+|P|²/m²)`. At a field point with lab coordinate `X=x-center`, the rest coordinates are

\[
t_0=\Gamma(t-v\cdot X),\qquad
X_0=X+{\Gamma-1\over v^2}v(v\cdot X)-\Gamma vt.
\]

The rest center moves with velocity +v in the lab and therefore has positive ADM momentum `m Gamma v`. These formulas specify the inverse coordinate map to be used for tensor pullback. Apply a proper spatial rotation to align the rest spin first; then apply the boost in the requested lab direction. First and second derivatives of the four-metric require all corresponding Jacobian factors. Setting `partial_t gamma=0` after a boost is incorrect even though the rest metric is stationary.

For a centered isolated unattenuated seed, asymptotic charges must satisfy `E=m Gamma`, `P=m Gamma v`. At the same coordinate origin, Lorentz transformation of the rest angular-momentum tensor gives `J_parallel=S_parallel`, `J_perp=Gamma S_perp`; these are coordinate ADM angular momenta, while the input S is the rest spin. Translating the center by c adds `c cross P` to J. A rotation Q transforms gamma/K as covariant rank-two tensors and transforms P/J as vectors. Only proper rotations (`det Q=+1`) preserve the spin convention directly.

The dissertation's 2015 parameter interface uses **lab ADM spin** as input and converts it to the rest spin before making the metric: printed pp. 131–133, Eqs. (5.1.21), (5.1.26), and (5.1.29) give `S_lab=Gamma S_rest-Gamma² v(v dot S_rest)/(1+Gamma)` and `S_rest=S_lab/Gamma+Gamma v(v dot S_lab)/(1+Gamma)`. The replacement API explicitly accepts rest spin. Explicit boosted thesis input S must therefore be converted before a like-for-like seed comparison; for perpendicular spin, divide the lab spin by Gamma. Tables reporting only dimensionless chi without explicit S/rest versus lab identification should not be silently converted or equated to rest S/m². A binary's measured horizon spin remains a separate observable.

A convenient axis-free QI Kerr four-metric uses `n=X0/r`, spin unit vector s, and `l=s cross X0`. Its spatial part is

\[
\gamma_{ij}={\Sigma\over r^2}\delta_{ij}
 + {a^2(1+2mR/\Sigma)\over r^4}\,l_i l_j,
\]

with `g_0i=-(2maR/(Sigma r²)) l_i`, `g_00=-(1-2mR/Sigma)`. This avoids spherical-axis divisions. It reduces to isotropic Schwarzschild at a=0.

## Conformal decomposition and operator identities

Physical reconstruction must use

\[
\gamma_{ij}=\psi^4 h_{ij},\qquad
K_{ij}=\psi^{-2}A_{ij}+\tfrac13\gamma_{ij}K,
\quad A_{ij}=M_{ij}+(\mathbb L_h b)_{ij}.
\]

All tildes in the reference correspond here to h. Covariant A scales with `psi^-2`; contravariant A scales with `psi^-10`. In particular the trace-free seed is

\[
A^*_{ij}=(\psi^*)^2\left(K^*_{ij}-\tfrac13\gamma^*_{ij}K^*\right).
\]

Any positive seed scalar that regularizes the puncture gives an admissible CTT decomposition, but different choices generally yield different *binary free data*. Generic Kerr in paper II.B.3 (PDF p. 7) chooses `psi*=det(gamma*)^(1/12)`, so each isolated conformal seed has unit determinant. The earlier analytical Kerr subsection instead uses `(Sigma/r²)^(1/4)`; boosted Schwarzschild subsection uses `1+m/(2r0)`. These choices must not be silently interchanged when comparing numerical binary values. The superposed h generally does **not** have unit determinant even if both seeds do.

Raise M using the *superposed* h inverse only after summation and trace removal:

\[
M_{ij}=M^{raw}_{ij}-{1\over3}h_{ij}h^{kl}M^{raw}_{kl},
\qquad M^{raw}_{ij}=A^+_{ij}+A^-_{ij}.
\]

Differentiating this projection includes derivatives of h and its inverse. Projecting each seed with its own metric and summing is insufficient.

For a contravariant vector b,

\[
(\mathbb L_h b)^{ij}=D^i b^j+D^j b^i-{2\over3}h^{ij}D_k b^k,
\quad \Delta_{L,h}b^i=D_j(\mathbb L_h b)^{ij}
 =D^2b^i+{1\over3}D^iD_jb^j+R^i{}_j b^j.
\]

The Ricci term is essential. Equivalently evaluate the divergence directly, with all derivatives of inverse h, connection and b retained. `D_j T^ij = partial_j T^ij + Gamma^i_jk T^kj + Gamma^j_jk T^ik`. The scalar operator is `h^ij(partial_i partial_j f-Gamma^k_ij partial_k f)`.

## Three different attenuation operations

The far filter is `F(r)=exp[-(r/s_far)^4]` at distance from the *same* seed. Paper Eqs. (17)–(20), PDF p. 7 apply it to `(h*_ij-delta_ij)`, `(psi*-1)`, and `K*`, while leaving covariant `A*_ij` unchanged. All derivatives must correspond to the filtered fields.

The companion filter is `f_±=1-exp[-(r_∓/omega_±)^p]`: each seed's metric and mean curvature are suppressed near the *other* puncture. Equations (22)–(25), PDF pp. 7–8 give

\[
h=\delta+f_+(h_+-\delta)+f_-(h_--\delta),\quad
K=f_+K_++f_-K_-,\quad
\Psi=\psi_++\psi_--1.
\]

Neither Ψ nor raw M is multiplied by f. The trace projection still follows metric superposition.

The equation filter g is the product of two C-infinity windows, with each window zero below r_min, one above r_max, and the tanh(tan(...)) interpolant in between. Test endpoints explicitly to avoid evaluating tan near its singular values. It alters equations only within its support. With `psi=Psi+u`, the later v3 Eq. (26), PDF p. 8 is

\[
F_H=\Delta_h u+g\left[\Delta_h\Psi-{\psi R_h\over8}
             -{\psi^5K^2\over12}+{A_{ij}A^{ij}\over8\psi^7}\right],
\]
\[
F_M^i=\Delta_{L,h}b^i+g\left[D_jM^{ij}
              -{2\over3}\psi^6h^{ij}\partial_jK\right].
\]

Do not take `D_j(g M^ij)` or multiply the correction operator by g: that changes the prescription. Do not multiply K by g *and* use Eq. (26) without documenting a new free-data choice. The earlier high-boost paper used a different prescription and also defined an explicitly approximate method which neglects some metric/trace terms; the latter leaves a nonzero exterior constraint floor and is not a substitute for the requested actual solve.

The dissertation confirms this historical distinction explicitly. Printed pp. 103–104, Eqs. (4.1.36)–(4.1.41), and pp. 135–136, Eqs. (5.1.38)–(5.1.44), instead put separate discontinuous step windows into raw `M=g_plus A_plus+g_minus A_minus` and `K=f_plus g_plus K_plus+f_minus g_minus K_minus`. A typical step radius is `.2 r_H`; the author reports that step windows gave the best independently measured charges among the tested windows. These equations also treat summed seed tensors as automatically trace-free and use the difference-of-connections divergence identity without the later full superposed-metric trace/index corrections. The v3 source-attenuation equations and explicit trace projection are the appropriate target for the new implementation; exact 2015 binary reproduction is a different formulation.

The optional operator modification (Eqs. (27)–(28)) uses `h_op=delta+g(h-delta)` and `Gamma_op=g Gamma` only in correction operators. This connection is intentionally not the Levi–Civita connection of h_op. Its extension in the transition region must be documented. Outside g<1, the actual h connection and physical vacuum constraints are required. Near g=0 the system becomes scalar/vector flat Poisson. Derivatives of the operator's chosen coefficient/connection extension must be treated consistently with the implemented operator.

## Infinity and puncture conditions

Physical boundary conditions are `u→0` and `b→0`, with finite regular correction fields at punctures. TwoPunctures enforces infinity by factoring `A-1` from each unknown. Do not impose `b~1/r²` universally: the balanced-momentum examples have this stronger falloff, but a net momentum correction can have `1/r` terms. Verify charges on successively larger spheres and fit inverse radius terms. Far-filtered backgrounds remove some leading logarithmic asymptotics but do not guarantee exponential convergence of all terms.

## Independent acceptance gates defined before expensive runs

These are proposed numerical gates for this replacement, not claims that the source prescribed these exact thresholds.

1. Analytic/manufactured derivatives, curvature and operators: scaled relative/absolute errors ≤1e-9 for regular double/extended-precision points; compare AD derivatives with an independently formed finite-difference sequence. Use nonunit determinant metrics. At a zero exact result use a fixed absolute scale, not a division by zero.
2. Single seeds: physical off-grid H and M errors below a validated independent FD floor, normally ≤1e-7 in mass-normalized units outside throats, with at least second-order/expected stencil convergence before roundoff. Exact Schwarzschild must have K=0 and vanishing three-Ricci scalar outside the puncture. Exact Kerr must have K trace≈0. Spinning/boosted seeds must pass positive-definite gamma, tensor symmetry, signed-normal continuity, and asymptotic charge errors ≤1e-3 after radius extrapolation.
3. Moderate binary solve: final nonlinear unscaled residual ≤1e-9 or documented stricter/tolerated precision floor; independently sampled exterior RMS of H and norm(M) ≤1e-6 and exterior maxima ≤1e-4, with systematic reduction across at least three spectral resolutions. A small internal residual alone is never sufficient. Record near-hole bins so a large bulk domain cannot dilute a local error.
4. Coordinate tests: relative tensor and scalar comparisons ≤1e-8 at corresponding rotated/translated sample points, and charge consistency within extraction uncertainty. These use different orientations and center offsets, not just permuted symmetric cases.
5. Publication-level high-spin/high-boost cases: require the moderate gates first. Aim for off-grid exterior RMS≤1e-6 near horizons and ≤1e-9 in a bulk shell after verifier error calibration; demonstrate convergence and report the observed limit if these goals fail. The paper's high-spin orbital near-horizon constraints reach roughly 1e-7–1e-6, and high-boost companion work has a near-axis H floor around 1e-7. It is inappropriate to demand source plots' bulk roundoff values from an uncalibrated lower-order verifier.

For every physical residual report both raw dimensionless `m_scale² H`, `m_scale² sqrt(gamma_ij M^i M^j)` and a normalized local residual with denominator `|R|+K²+KijK^ij` or the analogous term norms. Define any denominator floor. Preserve H and all M components, RMS/L-infinity, point coordinates, FD order and step, sample-set definition, resolution, solver tolerance, runtime, and iteration counts.

The independent verifier should consume only reconstructed raw physical gamma/K and differentiate them on off-collocation points. It should not reuse conformal R, seed derivative arrays, or the solver's CTT operator as its primary computation. For each sampled region refine verifier step `h,h/2,h/4`; increase solver resolution with verifier step fixed near its converged plateau. Otherwise derivative noise can masquerade as solver error or hide it.

Useful exact manufactured checks are `h=e^(4w) delta`, for which

\[
R_h=e^{-4w}[-8\Delta w-8|\nabla w|^2],\quad
\Delta_h f=e^{-4w}[\Delta f+2\nabla w\cdot\nabla f].
\]

For arbitrary smooth w, translations b=c, dilation b=x, and rotations b=omega cross x are conformal Killing vectors: Lb=0 and Delta_L b=0 despite nonzero connections. In flat space b=(x²,0,0) gives Delta_L b=(8/3,0,0). These checks expose both connection cancellation and nontrivial principal parts.

## Horizon enclosure and support

For an isolated QI Kerr seed `r_H=sqrt(m²-a²)/2`; boosting the seed contracts the minimum lab coordinate horizon radius to `r_H/Gamma`. A spherical r_max larger than this is certainly not inside the isolated seed horizon. This bound is a screening check, not a binary apparent-horizon certificate. At `m=.5, chi=.99`, the unboosted QI throat radius is only about .03527. Fixed source radii appropriate at moderate spin therefore cannot automatically be reused at near-extremal spin. LES/fisheye coordinates were introduced partly to avoid this small geometric scale.

Find the physical binary apparent horizon(s) and verify that the complete `g<1` support lies inside them before calling the exterior construction horizon-contained vacuum. Without a horizon finder, report attenuation-region and g=1-region constraints and state that horizon enclosure is unverified. A coordinate cutoff labelled exterior is not proof of enclosure. Horizon masses, irreducible masses and spins must be measured, not equated to seed m or S.

## Published initial-data targets and reproducibility limits

Paper Table I (PDF p. 10) provides useful unboosted high-spin parameters at centers `(±6,0,0)`:

| Case | Seed masses | Seed spins | Momentum | Reported ADM energy |
|---|---|---|---|---|
| HS90UU | .5,.5 | (0,0,.225) both | 0 | .982353 |
| HS90UD | .5,.5 | (0,0,±.225) | 0 | .982388 |
| HS99UU | .5,.5 | (0,0,.2475) both | 0 | .980124 |
| HS99UD | .5,.5 | (0,0,±.2475) | 0 | .980163 |

Fig. 2 uses omega=1,p=4 and large collocation counts. The exact inner/far filter choices and all numerical details are not supplied in Table I; these must be recorded as implementation choices, not inferred exact reproduction. Likewise Table III provides HSQC `b=4.7666,m=.48745, P_plus=(-.001138,.09794,0)`, opposite minus momentum, spins zero, ADM=.98914. Table V supplies unequal masses/generic momentum magnitudes and LES high-spin orbital parameters, but a QI replacement cannot expect exact matching to LES free data.

The high-boost appendix example specifies `d=66` and `P_x/m_H=±2`, but not a complete seed-mass/attenuation specification. Companion arXiv:1506.06153v2 Tables II–III supply measured momenta, horizon masses and separation for boosts through 4; these are physical targets, not complete seed files. They should not be described as exactly reproducible inputs until the missing free-data parameters are available. A fully specified local benchmark at a published dimensionless regime is useful, but must be called a corresponding benchmark rather than exact source reproduction.

The supplied dissertation does not unambiguously resolve the high-boost bare inputs. Printed p.119, Table4.3 gives `P/M_ADM`, measured `m_irr/M_ADM`, Gamma and separation. Its Gamma=`sqrt(5)` row is `P/M_ADM=.4530`, `m_irr/M_ADM=.2263`, `d/M_ADM=200`; the Gamma=`sqrt(17)` row is `.4908`, `.1217`, `400`. Horizon masses cannot simply be substituted for seed masses. Conditional reconstruction is possible: **if** table P is the input seed momentum and table Gamma is the coordinate boost of printed p.99 Eqs.(4.1.10)–(4.1.11), then `m=P/sqrt(Gamma²-1)`. This yields approximately `.226496` and `.122700` in ADM-normalized units for those two rows. It is a table-informed approximate benchmark with an explicit assumption and rounding uncertainty, not a unique historical parameter file.

There is a material ambiguity in that assumption: printed p.121 Eq.(4.5.46) defines Gamma for the radiation fit as `sqrt(1+(P/m_irr)²)`. The first five table Gamma values closely track the **rounded** P/m_irr column; e.g. `3.1717≈sqrt(1+3.01²)`, whereas a bare input ratio3 would give `sqrt(10)=3.16228`. The final row's Gamma=`4.1231` conflicts with P/m_irr=`4.03`, which would give `4.15222`; a source typo or mixed definition remains possible. Thus Table4.3 plus rounded Gamma/P/d cannot establish which mass was used in the boost map or fully specify the historical widths/stuffing. Original parameter files would resolve this. For unboosted high spin, thesis Fig.3.3 (printed p.85) used companion width `.2`, power4 and reports off-horizon RMS about `1e-8` by N≈112; later v3 Fig.2 uses width1 and much larger resolutions. This is another material free-data difference. Thesis Fig.4.1 (printed p.107) used widths1, power4, `P_y=.0848`, separation12, and reached physical RMS below `1e-7` above N≈80. Neither result supports publication-quality convergence at N≤28 without independent evidence.

## Apparent source typos checked by derivation

- The Cartesian Kerr A_xy expression on PDF p. 5 contains `cos(2theta)`; spherical-to-Cartesian transformation gives `cos(2phi)`.
- The stationary Kerr beta_phi formula in the QI section uses a radial symbol where the rational expression is in the Boyer–Lindquist radius R.
- Appendix B PDF p. 23 prints `E[4 uhat/r delta]=(1/(8pi)) integral uhat dOmega`. From its own ADM metric surface integral, the coefficient is `1/(2pi)`: constant uhat gives mass `2 uhat`, as Schwarzschild requires. Use the general metric ADM integral, which also handles angular background terms.

## Implementation audit

Initial source inspection (before solver/tests): `HiSpID_geometry.cpp`, `HiSpID_jets.hpp`, and `HiSpID.h` were reviewed. The axis-free rest metric, boost tensor pullback, signed boosted lapse, first derivatives of seed K using second-order jets, companion/far filters, superposed trace projection, and Eq. (26) residual/JVP coefficients agree with the derivation above. The optional operator is explicitly defined by lower-index `h_op` contractions of the chosen nonmetric-compatible connection, followed by raising and divergence with `h_op^-1`; this is an interior extension with correct exterior limit.

One numerical issue was reported immediately: the first implementation clamped the tanh(tan) inner window below z=.025 and above z=.975. At z=.025 the true window is approximately `9.19456e-12` and its z derivative `9.38476e-9`, so this created a small derivative/value jump rather than a clamp already below double precision. The parent repaired the cutoffs to `.01/.99`; the omitted values are then about `2e-28`, below double precision. This finding is resolved. A potential full-jet cache size of tens of GB at publication resolutions was also resolved by the compact approximately 2.4KB-per-point `Cached` representation.

The first seed routine rejected signed lapse magnitude below 1e-12 at the exact QI throat. This limitation has been removed by the regular rest-K/graph-slice formulation audited below; exact Schwarzschild and boosted Kerr throat samples now pass. Sampling at the unsmoothed puncture remains excluded, since background generation evaluates the full seed before applying g.

### Solver and physical-verifier audit

The completed `HiSpID_solver.cpp`, `python/hispid.py`, and `validation/physical.py` were reviewed read-only. Four-field point interleaving, coordinate chain rules, local/lab proper rotations, superposed-metric contraction of A, analytic JVP coefficients, right-preconditioned restarted GMRES, two-pass orthogonalization, true-residual rechecks, FD preconditioner transformation and line-search state restoration agree with the intended algebra. The optional interior operator has a documented literal nonmetric-compatible contraction. The independent verifier differentiates raw physical tensors; its Ricci derivative indexing, inverse derivative, mixed momentum tensor, and raising/norm conventions passed the checks below.

The solver stops on a **weighted** residual, with row multiplier `(sin(alpha) sin(beta))^6`. Its tolerance is not an unscaled or physical acceptance gate. Diagnostics preserve unscaled extrema, which must be assessed separately. The parent's `HISPID_STATUS.md` predeclared an **initial** moderate gate of physical RMS `<1e-4` and charge stability `<.5%`; this is weaker than the reviewer's proposed gate 3 above and is explicitly distinct from publication precision. The initial validation harness follows that RMS gate but omits a maximum/internal-physical gate and compares only last versus first resolution, which is not enough to establish systematic three-level convergence. These limitations were reported; there was no successful high-regime result at that point. High-spin samples starting at `.6 m` also do not probe the near-extremal horizon at approximately `.07 m`.

The reviewer ran only small seed/operator/reconstruction checks and read the parent's saved binary solutions; no independent binary solves were run. Checks used the explicit bundled Python with NumPy and one numerical thread.

| Independent check | Evidence |
|---|---|
| Raw Schwarzschild, Kerr, boosted Schwarzschild and generic boosted Kerr physical constraints at three off-grid points, including the inner sheet | On step `.004/.002/.001`, Schwarzschild H RMS `1.663e-7/1.037e-8/6.638e-10`; Kerr H `4.282e-7/2.681e-8/2.325e-9`, M `1.238e-7/7.737e-9/4.836e-10`; boosted Schwarzschild H `8.878e-8/5.582e-9/1.546e-9`, M `3.826e-9/2.396e-10/1.504e-11`; generic boosted Kerr M `9.595e-8/5.991e-9/3.742e-10`, H reaches a few `1e-9` roundoff floor |
| Signed boosted Schwarzschild mean curvature, including r below m/2 | Agreement with paper's analytic rational expression to `1.9e-15`; this checks K orientation, which vacuum constraints alone cannot establish |
| Generic isolated boost charges, radii 40/80/160 with inverse-radius quadratic extrapolation | E/P errors below `1.0e-6`, J error below `4.3e-6` relative to analytic Lorentz-transformed charges |
| Independent anisotropic manufactured physical metric `gamma=e^(4w) A`, `K=lambda gamma`, non-diagonal constant A | At step `.008`, H error `2.6e-11`, M error `6.1e-12`; M falls to `3.6e-14` by step `.002`. Expected `R=e^-4w[-8 tr(A^-1 Hess(w))-8 grad(w) A^-1 grad(w)]`, `H=R+6lambda²`, `M=-2e^-4w A^-1 grad(lambda)` |
| Independent curved vector identity `D²b+(1/3)grad(div b)+Ricci*b`, generic boosted Kerr conformal metric and arbitrary quadratic b | Versus native operator hook at `.008/.004/.002`: scalar error `9.8e-10/6.1e-11/3.8e-12`; vector error below `4e-10/2e-11/8e-11`; Ricci reaches FD roundoff around `1e-10` |
| Single unattenuated seed context, zero correction | Reconstructed gamma, K, psi, h, A agree with raw seed to `8.9e-16` |
| Constant auxiliary scalar U=.03, including map-axis samples | Sampled u agrees with `.03(A-1)` to `2.5e-17` |
| Random small JVP versus centered residual differences, 4×4×4 single-seed context | Relative difference `1.3e-8` at difference step `1e-3`; smaller steps show expected residual cancellation. This is preliminary evidence; the parent's larger-direction test supplies the stricter gate |
| Rotated/translated generic two-hole reconstruction with nonzero manufactured corrections | 6²×4 context, U scalar `.015`, lab vector `(-.02,.01,.025)`, rotation `.73` radians about `(1,2,3)`, translation `(.7,-.2,.4)`, including map-axis sample: relative gamma error `9.6e-16`, K/A `1.3e-13`, h `4.5e-16`; scalar/vector corrections agree to `2.1e-17`. This checks geometry/frame/interpolation covariance without a solve; converged binary covariance remains a separate gate |

### Confirmed sampling defect and failed binary gate

The saved moderate binary runs at 12²×8, 20²×12, and 28²×16 fail independent physical acceptance: near H RMS approximately `.029/.026/.059`, bulk H RMS `.005/.003/.009`; ADM energy extrapolations `-1.27/-.007/1.47` are unstable. Weighted residuals down to `1e-13` do not rescue this result. At the highest resolution the reported unscaled residual also rises as high as `.255`. High-regime acceptance must remain blocked pending repair and convergence.

The reviewer loaded the saved 12²×8 unknowns into a read-only context and checked six regular g=1 collocation points. Sampled corrections agree with stored `(A-1)U` to `1.8e-16`. The physical Hamiltonian relation `F_H=-psi^5 H/8` agrees with native residuals at the independent verifier floor: at the closest point, step `.004/.002/.001` gives differences `5.2e-7/3.3e-8/2.0e-9`. The large off-grid H therefore is real unresolved spectral structure rather than a value interpolation or Hamiltonian-contraction mismatch.

A concrete momentum-reconstruction defect was identified and reported: `make_coefficients` interpolates sampled phi derivatives from `Derivatives_AB3`. The inherited `fourder` correctly sets the cosine Nyquist mode's derivative to zero **at the collocation points**, but a trigonometric interpolant's derivative at arbitrary phi contains the missing term `-(N_phi/2) c_Nyquist sin(N_phi phi/2)/2`. Interpolating the derivative samples drops that term. Native `fourder2` retains the Nyquist second derivative. Consequently physical K constructed from the sampled b gradient is not the derivative of the sampled b interpolant, and its divergence disagrees even at collocation points.

For the saved 12²×8 solution, the closest collocation point had mapped physical momentum discrepancy `psi^10 M=(1.59e-5,-1.03e-4,-4.27e-5)` despite an approximately zero internal residual. The analytically predicted missing Nyquist principal term matched discrepancies at five other points to approximately `8e-15`–`3e-12`; the closest point differed only at its FD floor. The parent repaired the sampler by explicitly restoring the missing Nyquist derivative from its base coefficient. Independent recheck of the same saved solution confirms the persistent discrepancy is removed: the closest point's maximum momentum discrepancy decreases `5.96e-6/3.72e-7/2.33e-8` on steps `.004/.002/.001`, each ratio approximately 16; the other five points decrease `6.86e-10/4.29e-11/2.66e-12`. This finding is resolved.

A pure-Nyquist analytic test independently confirms the repair away from collocation: with a single unattenuated rest Schwarzschild seed (h=delta), eight phi nodes and auxiliary vector `U_b=c cos(4phi)`, the expected vector is `b=c(A-1)cos(4phi)` and `Kij=psi^-2[ci grad_j F+cj grad_i F-(2/3)delta_ij c dot grad(F)]`. At four arbitrary off-grid points sampled b and K agree to `7.5e-17` and `1.7e-16`. The repaired sampler also uses a rationalized prolate inverse and a four-transverse-point axis limit before inversion. The latter now agrees with the dissertation's printed p. 69 Eq. (2.4.212), including its `1e-4 M` default displacement; axis averaging is a finite-displacement approximation and should remain identified as such.

A useful scalar diagnostic for far filtering is an unboosted Schwarzschild configuration with h=delta, K=A=0 and g=1. The exact correction is `u=sum[(1-F) m/(2r)]`, restoring Brill–Lindquist psi and total mass. It independently tests the far-source shell, infinity factor, interpolation and scalar solve. The shell near F's scale 40 may need substantially higher resolution than 28; convergence must be measured, not assumed.

This diagnostic was run without solving, sampling the analytic correction into the native collocation unknowns for masses `.6/.4`, centers `±3`, far radius 40. At N=12/20/28, maximum psi interpolation error at seven fixed regular near/shell/far points was approximately `1.5e-3/5.1e-4/2.1e-4`; weighted scalar residual was `6.0e-4/2.2e-4/3.0e-4`. At N=64/96 with axisymmetric phi count 4, psi error fell to `8.2e-6/7.0e-7`, weighted residual to `1.3e-5/1.7e-6`. This confirms substantial underresolution of the far shell in the initial moderate grids. The diagnostic's global raw conformal residual is approximately `1e7` at the smallest coordinate radii, because small interpolation errors are amplified by the singular map.

That analytic control also clarifies the proposed internal gate: global absolute raw conformal residual `<=1e-9` is unsuitable at arbitrarily small puncture coordinate radii. Preserve raw extrema for diagnosis, but assess the internal residual with physical-equivalent factors `-8 F_H/psi^5` and `F_M/psi^10` in the actual-metric g=1 region, and report separate support/core and exterior bins. This is an explicit, analytically justified refinement of gate 3; it does **not** relax its independent off-grid exterior physical thresholds, maxima, or resolution convergence. A seed-horizon screening radius is still not a binary apparent-horizon certificate.

As a concrete check, the reviewer loaded the saved far-radius-zero 24²×12 moderate solution and selected 64 representative g=1 collocation points with coordinate radius below 20 and outside 1.1 times each isolated seed horizon radius. The physical-equivalent internal maxima were H=`2.36e-8` and M components `(5.52e-10,8.73e-10,1.90e-9)`, despite the global raw residual reaching `.145`. This confirms the importance of the distinction and shows that the weighted `1e-12` stopping tolerance does not automatically imply a physical-equivalent `1e-9` internal gate. These are collocation checks; independent off-grid errors remain the stronger convergence test.

### Exact far-filter scalar reparameterization and revised derivative backend

To avoid forcing the regular spectral unknown to resolve the far scalar cutoff shell, the parent proposed `u=W+v`, where `W=sum[(1-F)(psi_seed-1)]` is evaluated analytically with jets and only `v=(A-1)V` is expanded. This is an exact change of unknown, not a new free-data choice: keep filtered Psi, h, K and raw M unchanged, and reconstruct `psi=Psi_filtered+W+v=Psi_unfiltered+v`. W is O(r³) with vanishing first/second derivatives at a puncture, and O(1/r) at infinity, so v retains the admissible finite C2 puncture correction and infinity conditions. Generic `r³ C(direction)` is not necessarily C-infinity or analytic at the puncture, and v need not retain any stronger smoothness of u; this change alone does not prove exponential spectral convergence. The coefficient can be anisotropic for boosted seeds without invalidating the required first/second derivative properties.

Residual evaluation must include `Delta_op W` unmultiplied by g, including g<1, and all scalar W value/gradient/Hessian terms in the base field. At g=1, `Delta_h W+Delta_h Psi_filtered=Delta_h Psi_unfiltered`; the far shell cancels exactly. Newton/JVP directional fields do not contain W, but current-Newton/base fields do, including the public JVP wrapper. Sampling must report total u. The debug unknown scalar now represents V rather than the prior U; saved unknowns need parameterization/library metadata. Angular 1/r terms from boosted/determinant seeds may require cancellation by v at infinity, so this change is not itself evidence of exponential convergence.

The implementation was subsequently reviewed and independently checked. `fields(..., include_reference=true)` adds all ten scalar W jet entries to base fields; `jvp` passes false for directional fields. Physical sampling adds W to the reported u and total psi. A 6²×4 single-Schwarzschild test with far radius 2 and g transition `.1`–`.4` supplies an especially discriminating interior check: V=0 gives `F_H=(1-g)Delta W`, where `Delta W=(m/2)F[12r/s^4-16r^5/s^8]`. The returned mapped scalar equation agrees with this independent analytic formula to `9.2e-16`, with 24 g<1 points and expected maximum mapped H `.00125`. This verifies W is not incorrectly multiplied by g. Random JVP differences in this same context agree to `4.3e-11/3.9e-12/5.4e-13` on difference steps `.001/.0003/.0001`.

The sampler now differentiates its Chebyshev/Fourier coefficient basis directly in one traversal, including Nyquist derivatives naturally. An independent combined basis test with A, B, AB, A², sin(phi), sin(2phi), and cosine Nyquist terms compares physical psi/K and vector b against analytic Cartesian chain rules: errors `2.2e-16/4.9e-17/3.9e-17`. This checks the new implementation beyond the original single-mode repair.

`HiSpID_spectral.hpp` replaces repeated trigonometric transforms by cached differentiation matrices, without changing the collocation interpolant. Barycentric Chebyshev D and D², separate even-grid Fourier D2, stride/indexing, and mixed derivative actions were reviewed. The parent's native random-array test compares all nine derivative arrays to inherited `Derivatives_AB3` with normalized tolerance `1e-11`. The separate Fourier D2 is essential: squaring D1 would lose Nyquist second derivatives. The physical-equivalent equation hook rotates momentum into the lab frame and exposes g/psi; its header correctly states that mapped equations equal physical constraints only at g=1.

### Updated moderate and extreme-seed evidence

The parent's no-far-filter generic unequal-mass moderate run uses the same seed/companion/inner data at 24²×12, 40²×20, and 56²×28. Near H RMS decreases `3.25e-4/1.20e-4/1.50e-5`, and momentum norm RMS `2.62e-4/7.13e-5/9.12e-6`; at 56²×28 bulk H/M RMS are `1.18e-8/6.77e-8` and ADM energy extrapolates to `.9944889345`. These saved physical results were read, and pass the parent's declared **preliminary** moderate gate with three improving resolutions and stable charges. They do not yet meet the reviewer's proposed `1e-6` near-hole gate or publication accuracy. The original N56 solve took approximately 496 seconds before derivative-matrix optimization. A rerun with the regular graph seed and solved covariance checks remain pending.

The parent reports that the revised far-radius-40 series still fails systematic momentum convergence after the exact scalar reparameterization: near M RMS increases from `9.22e-6` at N56 to `1.45e-5` at N72, although N72 near H reaches `6.6e-6`, bulk H approximately `1.8e-7`, and ADM E `.99448074`. This failed series is retained. Passing the no-far-filter construction does not validate the filtered binary construction or justify high-regime filtered results.

The no-far-filter moderate and solved covariance reruns have since completed with the regular graph seed. Per-record library SHA matches the reviewer-tested library (`0016203beecd...`). Finest near H/M RMS remain `1.4973e-5/9.1217e-6`, with maxima `4.7792e-5/3.9933e-5`; bulk H/M RMS are `1.1380e-8/6.7693e-8`. The independent step-halving changes near H RMS by less than `.1%` and M RMS by approximately `.003%`, so these near residuals are spectral errors rather than verifier truncation. At g=1, mapped collocation maximum H is `2.35e-10`, all M components below `3.7e-11`, whereas modified-region physical H RMS is `.00137`. The finest solve took approximately 198 seconds with cached derivative matrices. Preliminary acceptance passes; the strict near-hole gate remains false.

Saved solved rotation/translation evidence was independently read. At N56, relative maxima for gamma/K/corrections are `4.1e-15/1.4e-13/3.3e-13`; translation tensor/scalar errors are at approximately `1e-13` or below, and the translated solve needs no Newton updates. The charge covariance gate passes its declared mass-normalized absolute `.005` threshold. Maximum charge/origin-fit differences plateau at `1.29e-4` (N40) and `1.25e-4` (N56), despite accurate field rotation; this is a finite-radius/quadrature accuracy limitation that should be reported and checked with increased angular quadrature before claiming precise charge covariance. The saved charge integrals use only 12×24 angular quadrature and quadratic fits at radii 100/200/400.

An additional small raw-seed check covers combined rest spin chi=.99 along `(1,2,3)` and lab velocity `(.86,.2,.1)` (Gamma approximately 2.2), without a binary solve. Near-throat and inner-sheet rest radii `.0687/.0291` show H errors falling from `(1.28e-7,1.82e-6)` at step `2e-4` to `(8.7e-9,1.14e-7)` at `1e-4` and `(6.5e-9,7.9e-9)` at `5e-5`; deepest-point momentum norm falls `8.50e-7/5.32e-8/3.31e-9`, with near-throat momentum reaching approximately `2e-8` FD/roundoff floor. Metric eigenvalues are positive. This extends evidence for seed geometry, **not** the validated binary parameter range.

### Regular graph-based extrinsic curvature (implementation audited)

The implementation removes the QI throat division by constructing extrinsic curvature on the boosted rest-coordinate graph `t0=-v dot X0`. The reviewer independently checked the derivation and prototyped it with finite-difference rest-coordinate derivatives, without changing source. Let beta=omega l with l the rotational Killing field; the stationary rest curvature is `K0ij=(l_cov_j diomega/alpha0+l_cov_i djomega/alpha0)/2`. Since `R_r=(r-c/r)/r` and `Delta=(r-c/r)^2`, the displayed analytical expression for `diomega/alpha0` has no throat division and remains finite at `r=sqrt(c)`.

For graph tangents `E_a=T_a^b e_b-alpha0 v_a n0`, `T=I-v tensor beta`, orthogonality gives `u=-alpha0 gamma^-1 v/q`, `q=1-v dot beta`, and graph normal `n=W(n0+u)` with `W=q/sqrt(q²-alpha0² v gamma^-1 v)`. The Gauss identity `nabla_e e=Gamma e-K n`, normal derivative `nabla_e n=-K e`, and torsion identity `[n,e_c]=a_c n+(dc beta/alpha0)e` reproduce the parent's Q/N expressions and `Kgraph_ab=-W[T_a gamma Q_b+alpha0 v_a N_b]`, with the stated signs. The final spatial boost B pulls both tensor indices back to lab coordinates. Signed alpha0 must be retained even though no inverse alpha occurs.

Prototype evidence: for m=1, S=`(.2,.3,.4)`, v=`(.4,-.2,.15)`, rest radii `1.4/.245/.539`, the regular rest K agrees with the existing metric-derived seed to `3.6e-15`. Graph K using independent fourth-order FD derivatives agrees at the inner-sheet point to `4.4e-10/2.8e-11/2.1e-12` on steps `.001/.0005/.00025`; outer points reach approximately `2e-13` floor. Boosted Schwarzschild graph K agrees to `2.3e-12`. This verifies signs/factors and the inner-sheet convention.

The graph derivatives must be taken in **rest** coordinates. The original seed x jets carry lab derivative slots after boost; using those slots directly as `db u`, `dc beta` or rest Christoffels would be wrong. Construct independent rest-coordinate jets and pull their derivatives back with B after graph algebra, or explicitly undo the derivative map first. Only K first derivatives are needed by the constraints; physical gamma still requires exact second spatial derivatives.

The updated source follows this order: independent rest-coordinate jets form gamma0, alpha0, beta0, K0, rest connection and graph Q/N; tensor indices are pulled back with the spatial boost; conformal scalar/trace/projection are formed; all jet derivative slots are then pulled back with the spacetime boost. The Q/N contractions and regular `gradomega` formula agree with the derivation. Differentiated jets are valid through first order for K, as required; their K Hessians are not a supported third-metric-derivative interface.

An independent exact-throat test reads the internal Seed jets through a read-only ABI binding, first confirming every tested value equals the public seed API. For m=1, S=`(.2,.3,.4)`, v=`(.2,.3,-.1)`, center=`(.13,-.2,.06)`, rest direction proportional to `(.3,-.5,.8)` and rest radius `r_H=.4213074886588179`, gamma has positive eigenvalues `(19.43,20.54,24.94)` and K is finite/symmetric. Fourth-order differences of public K values versus native K first derivatives give absolute errors `1.47e-6/9.19e-8/5.75e-9/3.54e-10/3.00e-11` at steps `.004/.002/.001/.0005/.00025`; the final scaled error is `1.93e-12`. Mean-K derivative error reaches approximately `6e-13` floor. Physical gamma derivative error similarly falls to `4.56e-11`.

Independent raw-tensor constraints at the exact throat and rest radii `.98/1.02` times r_H show fourth-order momentum convergence: exact-throat M norm `3.09e-8/1.93e-9/1.21e-10/7.41e-12` on steps `.004/.002/.001/.0005`. Exact-throat H falls `6.33e-8/4.06e-9/6.77e-10` before second-derivative cancellation dominates at approximately `1e-9`; neighboring inner/outer points show the same behavior. This removes the earlier exact-throat API limitation. It does not establish horizon enclosure for attenuated binary data.

A second exact-throat control combines rest chi=.99 along `(1,2,3)` with v=`(.86,.2,.1)` (Gamma=2.1801), m=1 and rest direction proportional to `(.3,-.5,.8)`. At `r_H=.07053368`, gamma remains positive and K derivatives finite. On steps `4e-4/2e-4/1e-4/5e-5`, scaled K-first-derivative errors fall `2.94e-7/1.84e-8/1.15e-9/7.18e-11`; raw physical M norms fall `1.95e-7/1.22e-8/7.63e-10/4.71e-11`, and H falls to its approximately `2e-9` independent FD floor. This is seed evidence only and must not be presented as an extreme binary solve.

### High-regime interpretation and charge quadrature follow-up

The planned HS99UU sequence 48²×12, 80²×16 and 112²×16 uses the thesis Table 3.1 (printed p.83) seed values m=.5, centers ±6, Sz=.2475 and ADM target `.980124`, with Fig.3.3's (printed p.85) companion width `.2`, power4, QI scalar and no inner/far filters. The modern superposed-metric trace projection is explicitly retained instead of assuming the summed seed tensors are automatically trace-free. Consequently ADM agreement within the predeclared absolute `2e-4` target is a reference benchmark with an identified formulation departure; it cannot certify exact historical free-data reproduction or horizon masses/spins. The proposed stricter physical RMS `1e-6` and maxima `1e-4`, using horizon-scale near points and verifier-step calibration, are reasonable replacement acceptance thresholds, distinct from the thesis's approximately `1e-8` RMS at N≈112. Systematic three-level convergence remains required. A Gamma=sqrt(5) boosted Schwarzschild binary is a local benchmark unless all missing historical bare inputs/stuffing are supplied.

The charge API was reviewed: ADM E differentiates raw gamma with a fourth-order Cartesian stencil of step `1e-4 R`; P uses `(Kij-K gammaij)n_j`, and J uses `(x-center) cross Pflux`, with flat surface measure. Gauss–Legendre mu and midpoint uniform phi agree with the independent Python convention. The native K trace is prescribed background K; the independent implementation recomputes the trace from gamma/K. Angular quadrature and finite-radius fits remain distinct errors.

The following diagnostic is intentionally supplied as documentation only: the reviewer is permitted to write this review file, while source/scripts remain parent-owned. It restores a retained checkpoint without solving or writing files, validates exact configuration/grid/parameterization/library compatibility, compares native versus independent raw ADM, and isolates angular covariance with an **exactly rotated sample callback**. Thus covariance error here measures integration error rather than a second solved interpolant. Run later as the only numerical job; default radius200 isolates angular refinement. Use `--radii 100,200,400` to add asymptotic fixed-origin checks, and compare increasing radii separately when refining extraction.

```python
# Save as a parent-owned validation/check_charge_quadrature.py, then run with
# PYTHONPATH=python:validation:examples OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# python validation/check_charge_quadrature.py --library /absolute/libHiSpID.so
import argparse, ctypes as C, hashlib, json
from pathlib import Path
import numpy as np
from hispid import Backend, Config, Hole
from physical import charges, extrapolate

p = argparse.ArgumentParser()
p.add_argument('--library', required=True)
p.add_argument('--case', default='moderate_far0_stable')
p.add_argument('--resolution', type=int, default=40)
p.add_argument('--radii', default='200')
p.add_argument('--quadratures', default='12:24,20:40,32:64')
args = p.parse_args()
root = Path.cwd()  # Run from the isolated native worktree.
data = json.loads((root/'validation/results.json').read_text())
rec = next(r for r in data[args.case]['records']
           if r['resolution'][0] == args.resolution)
backend = Backend(args.library)
sha = hashlib.sha256(backend.path.read_bytes()).hexdigest()
if rec.get('library_sha256') != sha:
    raise ValueError('checkpoint/library mismatch; do not reinterpret old coefficients')
if rec.get('unknown_parameterization') != 'u=W+(A-1)V, W=sum((1-F)*(psi_seed-1))':
    raise ValueError('unsupported checkpoint parameterization')
cfg = backend.config()
for name, _ in Config._fields_:
    if name == 'memory_limit_mib' and name not in rec['config']:
        continue  # Resource-only trailing field absent from archived old ABI.
    value = rec['config'][name]
    if name == 'hole':
        for i, hole in enumerate(value):
            cfg.hole[i] = Hole(**hole)
    elif isinstance(getattr(cfg, name), C.Array):
        getattr(cfg, name)[:] = value
    else:
        setattr(cfg, name, value)
if list(cfg.n) != rec['resolution']:
    raise ValueError('checkpoint/configuration grid mismatch')
n, _, nphi = rec['resolution']
checkpoint = root/f'validation/raw/{args.case}_{n}_{nphi}.npz'
with np.load(checkpoint) as saved:
    unknowns = saved['unknowns'].copy()
radii = [float(r) for r in args.radii.split(',')]
quadratures = [tuple(map(int, q.split(':'))) for q in args.quadratures.split(',')]
axis = np.array([1., 2., 3.]); axis /= np.linalg.norm(axis)
W = np.array([[0, -axis[2], axis[1]],
              [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
Q = np.eye(3)+np.sin(.73)*W+(1-np.cos(.73))*(W@W)
offset = np.array([.7, -.2, .4])
def rotate_charge(q):
    return np.r_[q[0], Q@q[1:4], Q@q[4:7]]
output = dict(case=args.case, resolution=rec['resolution'],
              library_sha256=sha, radii=radii, records=[])
with backend.create(cfg) as solution:
    solution.set_unknowns(unknowns)  # Deliberately no solve().
    def rotated_sample(x):
        values = solution.sample((np.asarray(x)-offset)@Q)
        for key in ('gamma', 'Kij'):
            values[key] = np.einsum('ik,nkl,jl->nij', Q,
                values[key].reshape(-1, 3, 3), Q).reshape(-1, 9)
        return values
    previous = None
    for nt, np_ in quadratures:
        native = np.array([solution.charges(r, ntheta=nt, nphi=np_)
                           for r in radii])
        raw = np.array([charges(solution.sample, r, ntheta=nt, nphi=np_)
                        for r in radii])
        rotated = np.array([charges(rotated_sample, r, center=offset,
                                    ntheta=nt, nphi=np_) for r in radii])
        expected = np.array([rotate_charge(q) for q in raw])
        item = dict(ntheta=nt, nphi=np_, native_EPJ=native.tolist(),
                    independent_EPJ=raw.tolist(),
                    native_vs_independent_linf=float(np.max(abs(native-raw))),
                    centered_rotation_linf=float(np.max(abs(rotated-expected))))
        if previous is not None:
            item['quadrature_change_EPJ_linf'] = float(np.max(abs(raw-previous)))
        previous = raw
        if len(radii) >= 3:
            global_origin = np.array([charges(rotated_sample, r,
                ntheta=nt, nphi=np_) for r in radii])
            fit = extrapolate(radii, raw)
            expected_global = rotate_charge(fit)
            expected_global[4:7] += np.cross(offset, expected_global[1:4])
            item['independent_extrapolated_EPJ'] = fit.tolist()
            item['fixed_origin_extrapolation_linf'] = float(np.max(abs(
                extrapolate(radii, global_origin)-expected_global)))
        output['records'].append(item)
print(json.dumps(output, indent=2))
```

### High-spin convergence failure and stable free-source divergence (implementation audited)

The parent reports HS99UU 48²×12,80²×16,112²×16 near M RMS `1.86e-5/2.03e-6/3.14e-6` and bulk M RMS `1.22e-8/3.05e-8/1.33e-7`, while near H RMS improves `5.27e-4/9.63e-5/7.27e-6`. These saved records were independently read. Both finest near H and M miss the predeclared strict RMS `1e-6`, and momentum lacks systematic convergence; neither published reproduction nor high-spin binary acceptance is established. Finest maxima H/M=`3.41e-5/1.36e-5` pass the maximum gate, illustrating why maxima, RMS and convergence must remain separate checks. ADM E converges to `.9802286764`, differing from the thesis target `.980124` by `1.0468e-4`: this passes the predeclared absolute energy comparison but cannot override the failed physical gates. Raw conformal core M extrema of order1–5 suggest cancellation to investigate, but do not by themselves establish the cause of off-grid error, because physical mapping suppresses them by psi^-10. An angular refinement control at N80 with Nphi24 was run before changing the source algebra.

For a seed conformal metric hs and covariant trace-free seed A, define stable differences directly from the filters:

```
dh = (f_self F_self-1)(hs-delta) + sum_other f F(h_other-delta)
di = h^-1-hs^-1 = -h^-1 dh hs^-1
dAup^ij = di^ik h^jl A_kl + hs^ik di^jl A_kl
dGamma^i_jk = .5 h^il(Ds_j dh_kl+Ds_k dh_jl-Ds_l dh_jk)
```

Here every inverse superscript in the products denotes the inverse metric, and Ds is the seed Levi-Civita derivative. The inverse difference follows by multiplying `h-hs` between the inverses, and the Aup difference is a telescoping tensor product, so both are exact. Metric compatibility Ds hs=0 gives the displayed exact connection difference. Consequently

```
Div_h Aup_h = (2/3) psi_seed^6 hs^-1 grad K_seed
             + Ds_j dAup^ij
             + dGamma^i_jk Aup_h^kj + dGamma^j_jk Aup_h^ik
T = sum_seed di:A_seed                  # trace_seed A_seed=0 exactly
Div_h M = sum_seed Div_h Aup_h - (1/3) h^-1 grad T
```

The first term uses the exact vacuum seed momentum constraint in either conformal choice. It is a justified analytic cancellation, independently supported by the raw seed constraint checks above. The projected-divergence term uses actual h's metric compatibility. These identities must be applied to the actual metric/source, not the nonmetric-compatible optional interior operator connection.

The required AD quantities are only dAup first derivatives, trace T first derivatives, and the seed connection value; metric second derivatives suffice for differentiated metric/connection jets. K/A second derivatives are not required and their differentiated-jet Hessians must not be interpreted as accurate third-metric derivatives. Use the same stable T in projected M values/first derivatives and in its source divergence. Unboosted stationary Kerr has K=0 exactly because beta=omega*l with l Killing and l·grad omega=0; set its K jet to zero **before** building seed A, with exact v²=0 rather than a small-velocity threshold. Returning rest K0 directly also avoids a redundant gamma/gamma^-1 contraction. These changes preserve the continuum free data; isolated, small two-seed, derivative/JVP, and independent physical controls are required before accepting numerical improvements.

The source implementation was subsequently audited. `seed_sum_source` forms dh from f/F directly, di with the exact inverse identity, du with the telescoping inverse product, covariant Dh with both seed lower-index connection terms, and the two connection-difference divergence contractions with the correct indices. It uses the unfiltered seed psi/K gradient in the vacuum identity and subtracts the actual-metric trace gradient. Stable T is used in both projected M and divergence; unboosted K is set exactly zero before A formation. The optional operator connection is not substituted into the source.

An independent read-only check compares the stable source to physical momentum computed only from sampled raw gamma/K, using zero spectral unknowns in a small 4³ context. For generic unequal-mass boosted seeds, both conformal choices and far radii0/8, it covers ordinary, near-hole and deep inner-sheet points. Native h values match the public sample. The mapped stable source `psi^-10[DivM-(2/3)psi^6 h^-1 gradK]` agrees with the independent physical momentum at approximately `1e-14` at the ordinary point; near-hole differences fall `1.27e-9/7.93e-11/5.01e-12` and deep-point differences `9.57e-9/6.00e-10/3.75e-11` on steps `.002/.001/.0005`. This verifies the reconstruction, source contractions and required derivative order independently. No binary solve was run by the reviewer.

The old-source N80 angular control, increased from Nphi16 to24, improves bulk M RMS `3.05e-8` to `1.03e-8`, but near M changes only `2.03e-6` to `1.93e-6` and H is essentially unchanged. It does not pass the strict near-hole gate. The old-source library is preserved separately with SHA `0016203beecd...`; its checkpoints must not be silently loaded under the revised source library.

Stable-source moderate and solved covariance rechecks were independently read, now under SHA `79e96b4c151a...`, with source case `moderate_far0_stable`; their preliminary gates pass and covariance's finest K relative error is approximately `5.8e-14`. The native/adapter configuration adds a matching trailing `memory_limit_mib` field, default2048 and accepted range16–8192. High-spin runs explicitly set4096; the conservative memory estimate is checked against this recorded per-context budget before allocation. This resource change does not change the equations.

The stable-source HS99UU 80²×16,112²×16,144²×16 results were independently read. Near H RMS is `9.63e-5/7.27e-6/1.77e-6`, near M RMS `1.83e-6/1.30e-7/8.72e-8`, and bulk M RMS `2.16e-10/1.29e-12/1.04e-13`. Analytic isolated cancellation restores momentum convergence and removes most of the previous noise/symmetry floor. The strict gate still fails due to finest near H exceeding `1e-6`; the JSON's generic `passed=True` currently denotes the preliminary `1e-4` gate, while `passed_strict=False` records this distinction. High-regime acceptance/exit status must honor its predeclared strict physical and charge criteria. ADM convergence alone cannot complete this milestone.

### Axis-value finite-offset bias (reported, repair pending)

A source-only manufactured control identifies an actual limitation of averaging **values**, as well as derivatives, at the four finite transverse offsets on the map axis. Take Schwarzschild masses `.6/.4` at centers ±3, disable f/F/g, and set the vector correction b^i=0 (the map parameter b=3). In a 10×6×4 grid let `z=A+1`, eta=.01 and

```
u(A)=eta*(z²+z⁴/6+23z⁶/720-49z⁸/1440),  U=u/(A-1).
```

U is a degree7 polynomial exactly representable by this grid; u vanishes at infinity. Since z=2 tanh(X/2), u is even in X and gives a smooth Cartesian field at the interpuncture axis. At x=y=z=0, analytic u=0 and `Delta u=4 eta/b²`, so physical `H=-8(4eta/9)/(7/6)^5=-.016450288570238586`. The public axis sample instead gives `u=1.11114e-11`, agreeing with the finite-offset bias `eta*(eps/b)²` for eps=1e-4. Independent H errors grow `3.21e-6/1.29e-5/5.14e-5/2.06e-4` as verifier h decreases `.008/.004/.002/.001`, a clear h^-2 center-value bias. No solve was run.

This does not explain the off-axis high-spin errors, but prevents strict physical constraint tests centered on the map axis. The proposed repair is to retain the documented finite-offset derivative limit while evaluating the central value analytically at its exact A/B coordinates using the Fourier zero-mode (azimuthal average) limit. Axis-gradient approximation remains a separate accuracy limitation. The manufactured control should be retained as a regression before claiming repaired axis constraint accuracy.

The precise control uses n=[10,6,4], hole0=.6 at (3,0,0), hole1=.4 at (-3,0,0), both spin/velocity zero, conformal_choice=0, far_radius=0, omega=[0,0], inner_max=[0,0], inner_flatten=0. Evaluate A_i=-cos(pi*(i+.5)/10). Broadcast U_i over Fourier and B indices into scalar channel of a C-order array shaped (4,6,10,4), with the three vector channels zero. The independent analytic point is (0,0,0); psi=7/6 and momentum vanishes.

A read-only callback simulating only the central-value repair, leaving the public sampler otherwise intact, gives physical H absolute errors at h=.032/.016/.008/.004/.002/.001 of 5.04e-9/3.12e-10/1.26e-11/2.29e-10/1.20e-9/6.63e-9. Coarse steps show fourth-order convergence before the expected double-precision derivative floor. Proposed native regression bounds are |u_axis|<=1e-14; H error<=2e-8 at .032, <=2e-9 at .016/.008 and <=2e-8 at the finer steps, with coarse error reduction .032/.016>=8 until the noise floor. These are calibrated from the simulated value repair, not a completed native regression; the actual repair remains to be tested after the current high-spin solves.

### Validation harness audit during high-spin refinement

The stage exit status now uses `passed_strict` for highspin/highboost, as required. Two remaining fail-closed details were reported to the parent: the high-stage prerequisites still recognize only the old moderate/moderate_far0 labels rather than the current stable-source series and do not enforce current-library SHA on the prerequisite seed/moderate/covariance evidence; and highspin's energy_agreement is recorded but not yet included in the combined stage acceptance. The current stable-source rechecks and energy comparison have separately been inspected, so these are reproducibility/harness defects rather than a reinterpretation of the existing physical evidence.

The parent saved the proposed read-only charge quadrature script. Its default case still selects old-source moderate_far0, which safely fails the SHA guard with the current library; the intended default is moderate_far0_stable. The optional compatibility skip should be limited to the resource-only trailing memory_limit_mib field, rather than silently skipping any absent configuration field. Native-only quadrature and incremental output options do not run a solve or change the retained coefficients.

### HS99UU azimuthal refinement and verifier-floor qualification

The parent completed 160²×16 followed by 160²×24 under the unchanged stable-source library. At N160phi16 the near H/M RMS is 1.55090e-6/8.68914e-8; increasing phi to24 reduces it to 3.17793e-7/2.69769e-9. The latter near maxima are 1.79295e-6/9.51411e-9. It passes the separately declared local strict RMS/max thresholds. Its bulk H/M RMS is 1.96005e-9/2.79910e-15; extrapolated E=.9802286761166793 differs from the thesis target by 1.04676e-4, within the predeclared2e-4 comparison tolerance. These records were independently read.

The aggregate strict gate remains false: the last three bulk H RMS values 1.79708e-9/1.53984e-9/1.96005e-9 are not monotonically improving, although the other tested raw norms improve. Preserve that original result. The newest stencil study gives bulk H RMS 4.30533e-10/1.96005e-9/5.96023e-9 at h=.004/.002/.001; smaller steps increase second-derivative roundoff. Near H changes only3.16796e-7/3.17793e-7/3.19610e-7 and near M2.71413e-9/2.69769e-9/2.69854e-9. This strongly supports, but alone does not conclusively establish, a bulk verifier floor.

A scientifically useful supplemental assessment may identify the bulk result as below calibrated verifier resolution, while retaining the formal monotonic gate failure. Required evidence is an analytic vacuum calibration at the identical point sets, magnitudes and stencil sequence, preferably both exact Brill–Lindquist with the HS masses/centers and a native exact chi=.99 Kerr seed at its actual center. A wider h=.016/.008/.004/.002/.001 sequence should locate the plateau independently of the solved residual; do not select the smallest solved residual as an estimate of truth. Pointwise differences between fine-grid bulk H should lie inside the measured noise envelope. A new binary solve is unnecessary for these read-only checks. No floor-qualified acceptance has yet been established by the reviewer.

The present normalized_M diagnostic is uninformative for maximal K=0 data: its denominator uses the norms of the already-contracted total DivK and gradK. Then DivK equals M and gradK vanishes, so normalized_M tends to1 whenever the residual exceeds the1e-8 denominator floor, even if the absolute momentum error converges. The HS normalized_M near .95 is therefore not a separate physical failure. Use the raw momentum norms as the gates, and either state this limitation or normalize by norms of individual derivative/connection contributions before cancellation.

One lower-priority mapping edge was also reported. With one active hole and a coincident inactive companion center, HiSpID_create's arbitrary-separation fallback leaves the active puncture at local x=0 rather than either map focus ±b. Existing controls intentionally use distinct companion centers, and both active binary cases are unaffected. The solver configuration should reject this edge or choose an artificial inactive focus that preserves puncture alignment. The standalone seed API does not use this map.

### Pending API-only repair source audit and high-boost status

The parent added an exact Fourier-zero-mode value at the central axis coordinates while retaining the four-offset derivative limit. Read-only audit confirms the .5 Fourier zero-mode factor is needed in addition to both Chebyshev T0 normalizations; coefficient indexing is correct. For |x|<=b it uses A=-1 and B=tan(acos(x/b)/2-pi/4); outside it uses X=acosh(|x|/b), A=2tanh(X/2)-1 and B=-sign(x). These are the correct prolate inverse limits. The new manufactured physical-H regression uses the precise regular degree-seven V control and the independently calibrated finite-step bounds. Coincident map centers, including an inactive companion focus, are now rejected in source. Neither change has been rebuilt/tested yet: the parent preserves library79 for the numerical jobs and retained checkpoint SHA guards.

The rho<1e-10*b branch retains a tiny-band approximation: nonzero transverse coordinates in that band are snapped to the axis correction limit. For a generic smooth component with a transverse gradient its value error is O(rho*|grad_perp|), not necessarily quadratic. Exact-axis values are the actual repaired limit. The offset eps>=1e-8*b exits the near-axis threshold and prevents recursion. API documentation should retain the finite-offset derivative accuracy and this near-axis approximation; the manufactured regression does not prove every solved field's finite-grid axis regularity.

The verifier-calibration source uses exact analytic BL and both native chi=.99 Kerr seeds on the actual18bulk points, then the saved last three HS grids across h=.016/.008/.004/.002/.001 without a solve. Its BL callback initially omitted attenuation, a required physical.constraints input; this was reported before execution. Calibration results remain pending.

The direct local Gamma=sqrt5, far_radius0 benchmark has not converged. Its 48²×8/80²×8 records were independently read: both native statuses are1 with Krylov iteration-limit errors, near H RMS .004929/.001916 and M RMS .018123/.004764, far above either acceptance gate. These failed iterates are retained. A N112 direct attempt and planned velocity continuation are provisional numerical investigations, not evidence of high-boost binary capability.

### Additional search for a fully specified published high-boost case

The review rechecked v3's complete relevant case descriptions, the actual arXiv TeX source archive, thesis chapter4's head-on sections, and the primary companion paper. No complete published high-boost bare configuration was recovered from these accessible materials. This is a scoped availability statement, not a claim that original author files do not exist.

| Candidate | Explicitly supplied | Inputs still missing or ambiguous |
|---|---|---|
| v3 AppendixA2c, PDFpp21–22, TeX2038–2134 | Equal nonspinning head-on binary; d=66M; Px/mH=±2 | Bare seed masses and hence input velocities; case-specific f widths/p; g radii; explicit case F and optional-operator selection. Horizon-mass curves are measured outputs. |
| Thesis pp113–119 d66/Px/mH2 case and Table4.3 | Equal nonspinning head-on setup; the chapter4 construction fixes p4, no farF, actual-metric correction operators and historical discontinuous stuffing | Bare seed masses; head-on f widths; exact stuffing radius (lambda=.2 is stated as typical); table Gamma's fit-vs-input ambiguity. |
| Companion1506.06153v2 II.2/Fig2, P/mirr1,d/M10 | Nonspinning head-on geometry, QI scalar, p4, specified smooth g functional form | Bare masses/input velocities, f widths, numerical g radii. Increasing g width is mentioned without a full parameter set. |
| Companion TablesII/III | Measured normalized momenta/masses, Gamma and separation | No complete bare masses or attenuation widths/radii. TableII is standard data; TableIII is approximate data omitting full trace/index corrections and cannot be treated as the exact coupled target. |

The arXiv v3 source tarball was inspected in memory: it contains spinorboost.tex, its bibliography and figure PDFs, but no implementation or .par files. Relevant hidden TeX comments also do not supply the missing high-boost parameters. The accessible public momentum_constraint branch is the previously reviewed BY/matter implementation. No original high-boost HiSpID parameter file has been found in the accessible source materials.

PaperTableIII / thesisTable4.1 HSQC is a useful bare-input orbital example, m=.48745 and P=(-.001138,.09794,0), centers±4.7666, but Gamma is only about1.020, so it does not satisfy a high-boost milestone. V3Fig6's explicit omega1,p6 applies to another case, Py=.0848,d12; it must not be silently attached to the different HSQC table binary.

A local fully specified benchmark or a conditional table-informed reconstruction is scientifically useful, but cannot complete exact published high-boost reproduction until the missing parameters are supplied. Do not substitute measured mirr or mH for bare seed mass.

### Completed verifier-floor calibration and current harness safeguards

The calibration JSON was independently read. At the identical18bulk points with h=.002, exact BL/Kerr control H RMS is1.158–1.298e-9, while fine-grid differences have RMS1.444/1.462e-9 and maxima2.585/2.680e-9; the control maximum reaches3.429e-9. At h=.016, control RMS is1.936–2.341e-11 and solved fine grids2.200–2.837e-11; grid differences1.903/2.292e-11 show no resolved spectral trend. The larger h is selected by the independent vacuum calibration, not by minimizing a solved residual. This supports the explicitly supplemental conclusion that near constraints converge and the tested bulk H is below calibrated verifier resolution. The original all-norms monotonic strict gate remains false. This evidence is bound to library79 and the recorded sample set.

Previously reported harness safeguards are now resolved in source: prerequisite gate records must all match the current library SHA; stable-source moderate labels and covariance's source case are checked; accepted_high_regime combines the strict physical/convergence gate with the declared highspin energy comparison. Charge replay now defaults to the stable case, restricts missing-field compatibility to memory_limit_mib and rejects ambiguous radial-resolution selections without Nphi. BL calibration supplies attenuation=1, and the native numerical comparison helper explicitly rejects NaN/Inf.

The target Gamma=sqrt5 exact-seed controls were independently read and pass physical H/M and refined charge gates. The earlier ADM momentum error1.15e-4 did not improve under angular refinement to32×64; shifting the radial fits from40–320 to160–1280 fixes the dominant asymptotic truncation. Both radial and quadrature sequences are retained, and the failed charge-resolution snapshot is preserved separately. The initial assumption of angular underresolution alone was not supported by the refined data.

For spinless collinear boosted seeds, rest beta=0 implies q=1 and the graph spacelike denominator is >=1-v²>0. The QI conformal boost-direction eigenvalue is Gamma²(1-alpha0² psi_seed^-4 v²)>=1, and transverse eigenvalues are1. With bounded f/F/g the actual and blended operator metrics remain positive. Thus there is no source/sign singularity expected around v=.89. The nonmetric-compatible optional operator can still have large shell derivatives g'*(h-delta) without the corresponding metric-compatible connection. It was suggested as a specific conditioning mechanism, not proven to be the sole failure cause.

A separate inner_flatten=0, far0 actual-metric-operator continuation reaches the target v=2/sqrt5 at48²×8 with7Newton/1120GMRES iterations and weighted residual max1.62e-15. Its full changed configuration and checkpoints are explicitly recorded. This establishes nonlinear convergence of that distinct configuration, not its physical accuracy: off-grid grid-convergence and charge gates remain pending. Earlier optional-operator Newton/Krylov failures remain retained; they must not be reclassified as a seed sign defect or proof that a solution does not exist.

The pending API migration source was audited. Default checkpoint/library guards remain strict; a separate explicit migration records both old/new SHAs, compares unchanged off-axis physical fields and residual/JVP, and characterizes several solved axis points without a new solve or automatic axis acceptance. The actual rebuilt migration/native/BY checks remain pending.

### Revised requested range and isolated AthenaK import audit

The user subsequently selects seed spin chi=.95 and boost speed v=.885 as the target range, adds an isolated AthenaK initial-data/horizon integration, and asks for measured solver acceleration afterwards. This supersedes insisting on Gamma=sqrt5 as the required boost target; the previous successful seed controls and failed binary attempts remain useful evidence. Here chi denotes rest seed |S|/m² and v the lab coordinate input speed until horizon measurements are available. Separate spinning-unboosted and nonspinning-boosted binaries do not establish the simultaneous chi=.95,v=.885 point. A joint range claim needs an explicitly specified combined spinning/boosted binary, including orientations and attenuation, plus the same independent physical and charge gates.

For m=.5, isolated QI Kerr horizon radius at chi=.95 is .0780625, compared with .0352668 at chi=.99. At v=.885, Gamma is2.147807737 and the smallest boosted coordinate horizon radius is .0363452. An inner-g radius chosen from the *spinless* radius can exceed that screen for a combined case; scale the support using the Kerr horizon radius and subsequently verify the actual binary horizon. The isolated screen is a planning check, not a binary enclosure certificate. New runs should predeclare any noise-qualified convergence assessment before inspecting their results, retaining the original raw monotonic flags.

The parent creates a separate AthenaK worktree at `/Users/hz0693/research/lazarus/.hispid-worktrees/AthenaK`, branch codex/hispid-pgen from PR790 head22baa243970fa1880b2bbc48e88a590069d55e47. The public [PR790](https://github.com/IAS-Astrophysics/athenak/pull/790) identifies the requested HengruiZhu99/project/z4c_overhaul branch. The reviewer reads that worktree's two_puncture pgen, id_solve HDF5 importer, ADM conversion, FastFlow, interpolation, and initial-driver path without changing AthenaK files.

The ADM storage is physical covariant gamma_ij and K_ij, with symmetric component order [xx,xy,xz,yy,yz,zz]. HiSpID full row-major tensors map through [0,1,2,4,5,8]. The old TwoPunctures pgen multiplies its *conformal* metric by psi^4; do not repeat that multiplication for HiSpID_sample's already physical gamma. ADMToZ4c uses c=(det gamma)^(-1/3), gtilde=c gamma, Khat=trace_gamma K at Theta=0, and Atilde=c(K-gamma Khat/3); it has the same K sign as HiSpID. Initialize the entire Z4c state to zero before conversion so Theta, shifts and gauge auxiliaries are defined. Z4cToADM reconstructs gamma/K with roundoff-sized algebraic projection effects. Its psi4=(det gamma)^(1/3) is the determinant-normalized Z4c factor, which generally differs from HiSpID psi^4 when det(h) is not1. GaugePreCollapsedLapse uses that ADM psi4 and is a freely chosen gauge, not the signed seed lapse used to derive K.

The checked id_solve importer reads metric/extrin [6,block,z,y,x], validates coordinate monotonicity, dimensions, ghost coverage, finite interpolated values and metric positive definiteness, and uses five-point tensor-product Lagrange interpolation. Its optional r_fill replaces interior gamma/K with a smooth BL-like mixture and must be disabled for an unmodified HiSpID import, or recorded as a further modification whose support is enclosed. The generic polynomial/layout tests do not validate near-horizon Kerr interpolation accuracy. A direct spectral checkpoint loader preserves the same native solution at every destination cell and avoids an extra Cartesian source-grid interpolation. Keep checkpoint/configuration/parameterization/library provenance strict; no solve is needed during import. Native context construction currently retains geometry caches, so the loader's memory and startup time should be measured separately from sampling and conversion.

FastFlow by default is scheduled only after an actual evolution stage. Driver nlim=0 or tlim=0 therefore does not invoke it. A dedicated initial-data-only path must explicitly call metric derivatives, Find(0,0) and Write(0,0) after the ADM fields are ready, set start_time=stop_time=0 for each horizon, and record that no evolution cycle occurred. Explicit centers with use_puncture=-1 avoid requiring evolved trackers. The default start_time is effectively infinity and stop_time=-1, so merely enabling num_horizons is insufficient.

Several upstream FastFlow defects or limitations are reported to the parent before using its output for scientific claims:

- SurfaceIntegrals and rr_min use Kokkos RangePolicy(0,nangles-1); that end is exclusive, omitting the final angular node. SurfaceIntegrals also leaves rho at zero there, biasing both the integrals and spectral flow. For an exactly round sphere the omitted area fraction is w_GL,last/(4*ntheta): .0118463 at ntheta5, .000355122 at17 and .000113938 at25. The corresponding mirr relative biases are approximately-.00594082, -.000177577 and-.0000569706. These values are independently computed from the same Gauss-Legendre weights, without an AthenaK run.
- The reported hrms is ∫H²dA/A, a mean square, without the square root. Reported hmean is ∫HdA without division by area. Expansion units and acceptance must reflect this or repair the diagnostics.
- ah_found depends on the absolute difference of successive mirr values being below mass_tol; it does not require a small expansion norm. The hmean_tol guard tests an unnormalized signed integral and cannot exclude cancellation. Require an explicit expansion RMS/max gate and radius/area/spin refinement before treating this flag as a converged AH.
- MetricDerivatives fills dg only at active mesh cells. MetricInterp uses 2*ng Lagrange nodes that include dg ghosts near block faces, but no derivative-ghost computation or exchange occurs. Loading gamma/K ghosts alone does not repair that path. A valid mesh-based finder must fill those derivative ghosts consistently, or use a direct initial-time spectral geometry callback.
- Missing havepoint nodes silently contribute nothing. Require unique global coverage of every angular node, positive finite radial function and gradient norm, and fail closed if any node is uncovered. A positive mean radius is insufficient for a radial graph with negative angular radii.
- Spin is the integral against flat coordinate rotation fields, rather than an approximate-Killing-vector solve. For stationary axisymmetric Kerr aligned with its axis it is a useful analytic control; for distorted/nonaxisymmetric/boosted binary horizons it is a coordinate-spin estimate, and its derived Christodoulou mass must be labeled accordingly. The coordinate vectors are not generally tangent to a distorted surface.

The implemented expansion contraction itself is correct for the requested convention: H=div(s)+K_ij s^i s^j-K. Its adjacent comment's u^-3 factor on the K contraction is a typo; code uses u^-2 as required.

An independent source-only interpolation control quantifies the derivative-ghost defect. For a block x∈[-1,0], nx8, ng2 and a linear positive metric gamma_xx=1+.01x, the exact derivative is .01 at every cell and point. At x=-.01 the inline four-node stencil uses active cells6/7 and ghost cells8/9; if unfilled ghosts are zero, it returns .00586496. At x=-.1 it returns .010595; at x=-.3 its fully active stencil returns .01. This is a finite mesh-face error, not a small roundoff limitation. The calculation reproduces the source's cell locations, indices and Lagrange weights without an AthenaK execution.

For accurate initial horizons, the preferred comparison is direct checkpoint-to-ADM mesh sampling followed by a mesh-refinement study, plus a direct FastFlow initial-time callback sampling gamma/K and gamma first derivatives from the same native solution. The latter removes mesh interpolation/FD error from the finder and provides a separate test of the mesh import. HiSpID's current public Point interface has no gamma derivative field; derivatives can be added from existing jets or calibrated independent Cartesian sampling. Validate the finder first on exact Schwarzschild, stationary chi=.95 Kerr, and boosted exact single Kerr/Schwarzschild surfaces. At least two angular orders, flow tolerances and initial guesses should give stable results. Certify attenuation enclosure using the full converged radial shape and actual support ball about each hole, not only a seed-radius screen or a single angular minimum about a mismatched center. These integration checks remain pending; the reviewer has run no AthenaK or binary solve.

For an exact single Kerr control, area=8πm²(1+sqrt(1-chi²)), mirr=m*sqrt((1+sqrt(1-chi²))/2), and rest S=m²chi. The boosted exact horizon has r_rest=rh, hence an ellipsoidal lab radial graph, with the same area. The coordinate-rotation spin diagnostic need not equal the invariant rest spin on a boosted/distorted surface; use area, expansion and known shape as the robust boosted control. A sufficient enclosure test with differing centers is |c_AH-c_hole|+r_modified < min R_AH, evaluated with angular refinement and a nonzero margin.

An independent angular representation control uses the exact boosted radial function R(mu)/rh=[1+(Gamma²-1)mu²]^(-1/2), with boost aligned to the polar axis. Orthogonal Legendre projection using160Gauss-Legendre points and a10001point test grid at v=.885 gives maximum relative shape errors .08254/.01115/.001493/.0001993/3.538e-6/6.270e-8 at lmax4/8/12/16/24/32. Even the exact surface is poorly represented by upstream's default lmax4, and lmax12 has about.15% shape truncation. These are representation errors, not measured FastFlow errors. They justify angular refinement and adequate nonlinear quadrature in boosted single-hole controls before binary horizon claims.

### Replay metadata defect discovered during integration

The parent reports that moderate_far0_stable records serialize memory_limit_mib=0, apparently copied from legacy structure padding, while the current library rejects0. The charge replay correctly fails closed. The reported physical coefficients and parameters may remain unchanged, but the complete recorded configuration is not directly replayable and must not be called fully reproducible until corrected. Preserve the malformed historical JSON. An explicitly documented resource-only migration to a new default2048 configuration, checked against the same physical inputs/library and compared through residual/off-grid fields, can validate retained coefficients without treating the original configuration as valid. The parent instead plans renewed three-grid evidence from the bound checkpoint with a factory-generated usable configuration; those results are pending and should be distinguished from the archived record. Default guards must remain strict; no automatic0-to-default conversion is justified.

The parent also reports an actual-metric-operator Gamma=sqrt5 N112 run with near H/M RMS1.24e-5/1.35e-3. Even though nonlinear weighted convergence succeeds, independent physical momentum strongly fails the preliminary gate. It does not establish binary high-boost accuracy. The revised requested chi=.95/v=.885 and combined cases remain the target, with all previous failed outputs retained.

### Read-only preparation for solver performance work

No performance implementation or benchmark has yet been reviewed. The existing source has several clearly separable phases: cache creation, spectral differentiation, Cartesian field transformation, point equation evaluation, sparse preconditioner row assembly/factorization/triangular solve, Krylov vector operations, and sampling/coefficient preparation. Measure their time and peak memory on an already validated case before selecting an acceleration backend, with identical tolerances and retained physical checks.

SpectralDerivatives::along has independent output rows; field transformation, residual/JVP point evaluation and preconditioner row construction are also independent. They are plausible OpenMP/SIMD targets. The mixed derivative passes depend on completed first derivatives, and Sparse::factor/solve have triangular dependencies, so preserve these orders/barriers and leave triangular work serial initially. The lazy coefficient preparation in HiSpID_sample mutates coefficients_valid/coefficients: initialize once before parallel sampling or guard that mutation. Worker exceptions must be caught inside the parallel region and propagated to the caller; thread_local last_error alone does not convey a worker's failure to the calling thread. MPI replication would multiply large geometry/Krylov caches and needs a measured memory plan.

A portable AthenaK checkpoint format should serialize explicit configuration fields and array shape/count, not raw Config bytes/padding. Record the W parameterization, input-frame rule and source/library provenance. A Perlmutter build naturally has a different binary SHA; require a distinct documented cross-platform migration with field/operator comparisons rather than weakening the same-binary replay guard. A direct FastFlow MPI geometry callback must give each angular node unique rank ownership before the existing reductions.

### Sampling-only/gradient API source audit and bounded FastFlow proposal

The parent adds HiSpID_create_sampler and HiSpID_sample_with_derivatives. Read-only source audit confirms the sampler avoids collocation geometry/derivative/Newton allocations, preserves the same origin/frame mapping and coefficient basis, and rejects solve/residual/JVP/equation-sample operations. The gradient formula includes both bg.psi.first and W.first plus the sampled correction gradient; lab pullback has all three frame factors. The derivative layout is [point][derivative][i][j]. Tests include rotated/translated exact BL with a nonzero far filter, equality of ordinary/full-context and sampler fields, and a nonflat boosted geometry gradient compared with independent fourth-order value stencils. Runtime results remain pending.

The initial64B/point sampler budget estimate is not conservative: values and coefficients already use64B/point, and make_coefficients calls legacy SpecCoef, which allocates two complete scalar3D temporary arrays, their pointer tables and a line vector, with inclusive n+1 extents. Peak heap usage exceeds80B/point. This is reported for correction before calling the budget conservative; a scratch-aware bound or a larger128B/point allowance is appropriate.

The proposed optional initial-time host callback supplies physical gamma[6], K[6] and dgamma[18] in derivative-major symmetric order. Rank0 may own every angular point while other ranks have local havepoint=0; retain the existing integral/coefficient reductions. Check global ownership using a **separate** reduced count array: do not overwrite local havepoint with global1, which would make nonowners integrate empty geometry. Owner data should be finite/SPD and callback failures must be caught and synchronized before collectives, preventing root-only exceptions from leaving other ranks waiting. The callback must be reset before its native context is destroyed and before any evolution; initial ADM mesh data then remain the normal evolution source. Start/stop times are stored in FastFlow at construction, so configure them before creation or explicitly update the object.

Find/Write also have a stale-result edge: FastFlowLoop clears ah_found but leaves the previous finite ah_prop values on failure, so Write can print an old horizon after an unsuccessful search. Clear result properties at each actual attempted search, while preserving last_a0 for initial guesses, or require explicit fresh-result/found status in every consumer. An initial-time path must record iteration/time, actual callback use and found/expansion/coverage status. Reporting sqrt(mean-square expansion) separately can preserve the historical hrms column's semantics; documentation must identify both units.

Minimal independent integration controls are a known flat round sphere R4 (area64π, H=.5, mean-square=.25, trueRMS=.5, spin0), exact Schwarzschild from inside/outside initial guesses (rh=m/2, area16πm²), stationary chi=.95 Kerr from .8/1.2rh guesses (known invariant area and signed S_z=m²chi at angular-refined quadrature), and exact v=.885 Schwarzschild with angular orders16/24/32 (known ellipsoid and invariant area). They require no nonlinear HiSpID solve: use zero unknowns with f/F/g disabled and a distinct inactive map focus. Intentionally loose mass stabilization on the flat non-horizon surface is a useful demonstration that ah_found must be supplemented by expansion acceptance. The actual callback/controls remain to be implemented and run by the parent; reviewer work stays read-only apart from this document.

### Rebuilt API migration, solved-axis characterization, and loader implementation audit

The parent reports all15 sampler/gradient/native-adapter tests passing and increases the sampler allowance to128B/point. **The original same-process axis_api_migration comparison is INVALIDATED and its zero-difference compatibility conclusion is withdrawn.** Subsequent investigation proves macOS dlopen deduplicates the old/new archives by their shared libHiSpID.so install name: both paths returned the same handle/symbol address, so the comparison did not load two distinct implementations. The historical JSON is preserved as invalid evidence. Single-build native/Athena controls are separate; a fresh-process comparison with verified loaded paths is required before any cross-library acceptance claim.

Solved-axis accuracy is still limited despite the exact manufactured axis regression passing. HS160²×24 axis H RMS at h=.008/.004/.002/.001 is8.9097e-7/3.5597e-6/1.4230e-5/5.6950e-5, almost exactly h^-2. Its M remains about6–7e-10. The moderate solved axis H/M is worse and also fails the exterior accuracy thresholds. These records do not justify an axis physical-accuracy claim.

Every nonzero azimuthal Fourier mode of a smooth **Cartesian component** must vanish as rho^|m|, with the remaining factor smooth and even in rho; m0 must be even in rho. This applies at interpuncture A=-1 and the outer B=±1 axes away from the punctures. The current finite Chebyshev/Fourier space does not explicitly enforce these endpoint conditions. Constant angular endpoint leakage produces an O(delta/h²) transverse Hessian, consistent with the observed factor4 sequence; an odd-rho m0 term would instead give h^-1. The central Fourier zero mode is then a specified axis value, not a unique limit of the finite-grid interpolant along every direction. A point-only averaging tweak cannot repair continuity. A durable regularity-enforcing basis or tau/boundary formulation would require new solver validation. C² blending around an exterior axis cylinder changes the physical data and must not be silently called the same vacuum solution. Avoiding exact-axis quadrature nodes with even ntheta is useful for the present x-aligned map, but does not prove uniform near-axis accuracy or generic rotated-axis behavior.

The new AthenaK pgen/reader/exporter and FastFlow extension were read. The pgen uses physical tensor order [0,1,2,4,5,8], initializes Theta/gauge auxiliaries, checks metric positivity, imports ghosts and compares an independently allocated saved ADM mirror with the ADM/Z4c round trip. Native gamma gradients are correctly reordered into derivative-major six-component FastFlow storage. RAII clears geometry_source and initial_shape before the native context is released, including on exceptions. Rank0 callback data and synchronized failure flags feed the existing reductions; global ownership is checked in a separate count vector while local havepoint is preserved. Both half-open bounds are corrected. The initial shape is projected using orthonormal real harmonics (m>0 includes sqrt2), so no additional normalization factor is needed. Shape coefficients are written with17digits. Found now additionally requires the opt-in sqrt(mean-square expansion) threshold, and attempted searches clear stale result properties.

The portable text format uses explicit ordered fields,17digit doubles, a version/parameterization/source-SHA tag, exact count, END and extra-data rejection. The source-library SHA is evidence metadata explicitly named in the Athena input; it is not by itself a checksum of the consumer executable. The exporter initially assigns strong/preliminary from the entire case aggregate even if --resolution selects an unaccepted coarse record. Require the selected record's corresponding local/strict gate and successful diagnostics as well, otherwise export only as an explicitly allowed diagnostic.

Additional fail-closed issues are reported before runtime acceptance: hispid_initial_horizons with no finders currently succeeds without doing a search; map configured active-hole indices rather than assuming hole0 when it may be inactive; ensure each finder time window includes0, because Find can return before clearing old state; require finite positive area/min radius and finite small expansion in addition to ah_found, preventing a skipped stale found flag from passing. The round-trip max accumulation must explicitly reject each nonfinite value because std::max can ignore a NaN comparison. Initial-shape projection also needs lmax>=1 and ntheta>lmax to avoid flow division by zero and unresolved/Nyquist harmonics. These are reported source findings; the parent is building/running the controls, and the reviewer runs no AthenaK jobs.

### Resolved loader guards and independent stable-harmonic audit

The latest source resolves the previously reported exporter/loader guards: selected checkpoint resolution must pass its own acceptance gate; diagnostic import requires explicit authorization; empty finder lists, inactive-hole ordering, missing time-zero windows, invalid shape resolution, nonfinite ADM round-trip components and nonfinite/nonpositive horizon properties fail closed. The text reader checks the 128B/point sampler allocation bound before allocating/reading unknowns. Athena_SINGLE_PRECISION is explicitly rejected for this import validation. Callback rank ownership, synchronized failures, derivative ordering and RAII reset remain correct on read-only reinspection. The ordinary mesh derivative-ghost limitation and solved-axis regularity limitation remain unresolved and must remain in the report.

FastFlow now uses a normalized scalar associated-Legendre recurrence, with the Condon–Shortley phase, rather than the general factorial Wigner sum. The reviewer independently evaluates the proposed equations in a short read-only NumPy calculation, with no AthenaK execution. The diagonal recurrence starts P00=1/sqrt(4pi), Pmm=-sqrt((2m+1)/(2m)) sin(theta) P(m-1,m-1); the upward-l coefficients are correct for normalized complex spherical harmonics. The theta derivative is [l cos(theta) P_lm - sqrt((2l+1)(l²-m²)/(2l-1)) P(l-1,m)]/sin(theta), taking the lower term as zero for l=m. The second derivative follows the associated-Legendre ODE. Explicit low-l checks give Y10=sqrt(3/4pi)cos(theta) and Y11=-sqrt(3/8pi)sin(theta)exp(i phi). FastFlow's existing sqrt2 multiplication for m>0 therefore preserves orthonormal real cosine/sine modes and the InitialShape projection normalization.

Independent numerical controls use seven angles [.025,.17,.7,1.1,pi/2,2.3,3.11] and l=[0,1,2,8,16,24,32,48,64]. The addition theorem relative maximum is1.19e-13 through l64. For all m at each tested l, 98-point Gauss–Legendre integration of 2pi P_lm² differs from1 by at most2.64e-13. Fourth-order theta finite differences at h=1e-4, angles [.17,.7,1.1,1.6,2.3,2.9], and selected pairs (l,m) through(48,24) agree in the first derivative to9.7e-12 and second derivative to1.8e-8 after normalization by max(1,max derivative magnitude); the latter is a finite-difference cancellation floor. An independent m0 control uses numpy.polynomial.legendre.legval and legder, rather than the C++ recurrence: at l32 the maximum absolute errors in value/first/second derivative are1.28e-15/4.80e-14/1.22e-11. These controls audit values and both theta derivatives; phi/mixed derivatives are the analytic cosine/sine differentiation factors.

Reproducing the *old* double-precision factorial m0 expression at theta=[1.1,pi/2,2.0] gives maximum absolute errors2.82e-8 at l32,3.98e-6 at l40,.00530 at l48, and about320 at l64 against the normalized Legendre recurrence. This demonstrates a material high-order harmonic cancellation defect and supports the scoped FastFlow repair. The new theta formulas divide by sin(theta) and are tested for Gauss–Legendre interior nodes; they do not implement exact pole limits. The general spin-weighted-harmonic implementation is unchanged.

### Initial-time exact-geometry integration evidence read by reviewer

Saved controls-basic-fixed and controls-flat-fixed JSON and representative run logs are independently inspected. Under executable SHA1335f813d7541e49e7d06640dc4a206ad21f02d70cf96401010a178feab528d8, exact Schwarzschild from .8/1.2rh guesses passes with relative area errors1.27/1.33e-14 and expansion RMS9.36/9.71e-8. Stationary chi=.95 Kerr lmax8/12/16, ntheta16/24/32 and distinct guesses passes with area errors<=4.11e-14 and signed Sz error<=1.69e-14; the lmax16 Sz error is5.33e-15. A representative Kerr run records ADM/Z4c round-trip error5.73e-16, time=0, cycle=0 and MeshBlock-cycles=0. The flat R4 control passes its deliberately loose non-horizon acceptance with area error2.22e-15, mean-square expansion .25, RMS .5 and spin0; its import round-trip error is0 and the log likewise confirms zero evolution. Earlier flat-fixture hmean failures remain preserved. The mixed basic JSON aggregate remains false because it includes the older flat-fixture failure, even though its Schwarzschild/Kerr rows pass; do not silently replace that historical aggregate.

The stable-harmonic boost885 series is also read, under executable SHA98e0466eac0425de6c5ccd677b5458044faa694dcba1bafbba8cd69d0ca2c651. Exact nonspinning m1 seed, v=.885, lmax16/24/32/48 and ntheta18/26/34/50 yield expansion RMS.0038387/.000120457/.0000100119/9.99233e-8. Coarse searches use a deliberately diagnostic .1 expansion tolerance; their strong row flags remain false. The finest passes the1e-7 expansion gate, has area relative error1.55e-15 and sampled independent radial-shape relative maximum3.78e-8. All four logs are verified by the harness to have zero evolution cycles. The independent written-shape reconstruction uses Legendre-polynomial derivatives and explicit normalization/Condon–Shortley factors, with the exact boosted ellipsoid as reference. Its maximum is measured on the stated off-grid sample set, not a certified continuous supremum. Fixed-l quadrature refinement and the combined spinning/boosted exact control remain necessary before claiming the integration's full requested range. The direct callback does not establish mesh-based finder accuracy or repaired solved-axis constraints.

### Mesh-constraint diagnostic and remaining integration edge checks

The newly added initial ADMConstraints path is read-only audited against src/z4c/z4c_adm.cpp. I_CON_M is the squared physical norm gamma_ij M^i M^j, so sqrt(sum(dx1 dx2 dx3 M2)/sum(dx1 dx2 dx3)) and pointwise sqrt(M2) maxima are correct. Hamiltonian and momentum contractions use the physical ADM metric/K and match the native K sign; the optional matter terms are absent for the vacuum fixture. This diagnostic differentiates the imported physical fields, including their loaded ghosts, rather than FastFlow's partially filled dg array. Thus the known ordinary mesh-finder dg-ghost defect does not apply to this ADM diagnostic.

Three separate masks report g=1, g<1 and a stencil-safe outer set. The guard ng*sqrt(dx1²+dx2²+dx3²) is conservative for Cartesian mixed-derivative offsets. Keep min_radius fixed above the largest grid guard across a refinement series to compare a common physical region; report the actual opt.fd_stencil as well as ng, since upstream dispatch permits different stencil choices. The current weights are coordinate cell volume dx1 dx2 dx3, not proper volume sqrt(det gamma) dx1 dx2 dx3. Label that convention explicitly. The g=1 mask alone can include stencils crossing g<1, and its residual must therefore remain distinct from the guarded set. The exact seed/finder controls do not substitute for an imported solved-binary mesh refinement study.

Two remaining source edges are reported to the parent. First, a binary with two active holes and only one finder can currently complete hispid_initial_horizons after checking hole0; excess finders are rejected but insufficient finders are not. A complete component-horizon enclosure certificate requires one surface for each active hole, or a separately supported common-horizon certificate; a partial search must be labeled partial. Second, the control harness reuses run directories while FastFlow appends summary/shape files. The summary reader takes the last row, but shape_error initially flattens the complete shape file and reconstructs from the first coefficients. A fresh run directory or strict last-row/expected-count check is required to avoid pairing a new summary/log with an old shape. Current freshly created stable-boost outputs are unaffected. The source/library SHA label in a checkpoint is provenance metadata, not a checksum of the consumer native library; record the explicitly resolved linked library, include path and consumer SHA, especially after changing HISPID_ROOT in a CMake cache or migrating platforms.

The latest source now resolves both integration edge checks above: component finder count must equal active-hole count, and the written-shape reader uses the last row with exact (lmax+1)² coefficient count. Mesh diagnostics explicitly record coordinate_cell_volume and the actual fd_stencil.

### Combined exact Kerr/boost flow stability control

The new combined chi=.95,v=.885 exact single-hole control at alpha1 initially fails at lmax16/24 after600iterations. This failure remains separate from the successful spin-only and boost-only controls. Read-only inspection of l24 verbose output shows initial area32.982401, already close to the analytic32.980437, and minimum sampled radius.0742464, close to the exact contracted radius.0726904. The unnormalized hmean envelope falls .140→.050 over about18iterations, then develops rapidly growing alternating signs and settles into a two-cycle (late areas34.0970/34.5690 and radii.03945/.02818). The final expansion RMS is about4.54. This is evidence for explicit-flow instability; it is not by itself evidence against the seed geometry.

The expected exact Kerr QI horizon is r_rest=rh=.5m sqrt(1-chi²), independent of rest polar angle. Stationarity permits the lab time-zero graph's varying rest t. The spatial boost map X_rest=B3(x-center) yields R(n)=rh/sqrt(1+(Gamma²-1)(vhat·n)²), independent of spin orientation, while the induced area is the invariant Kerr area. The horizon expansion, invariant area and known shape are the appropriate controls; the coordinate-rotation spin estimate is not required to equal rest spin on this boosted cut.

UpdateFlowSpectralComponents sets beta=.5alpha, B=.5 and A=alpha*(.5+1/[lmax(lmax+1)]). Lowering flow_alpha_beta_const therefore scales every explicit update linearly without changing the elliptic preconditioner. The flow's sigma normalizes the trace of the angular principal tensor: locally the corresponding angular eigenvalues are those of2G_AB/tr(G), positive and summing to2. An isotropic sphere gives1,1; a strongly anisotropic boosted Kerr surface can approach0,2, putting alpha1 near the alternating-sign high-mode stability boundary, with lower-order terms, mode coupling or nonlinear quadrature alias able to destabilize it. Throat size alone does not establish the cause. Alpha .1/.2 and fixed-l quadrature refinement are controlled ways to distinguish step stiffness from geometry error. At reduced alpha, a stable slow mode may need more than the original600iterations to reach1e-7; an iteration-limited monotone attempt must remain distinct from a divergent one. The parent owns these numerical controls; the reviewer runs none.

The preserved alpha.2 rerun subsequently passes the combined exact-seed control. The independently read controls-combined-alpha02/controls.json (consumer executable SHA9b3b059f116b3d9982b963e297d30cfc10faaf6af20b98e667b0d0c7a79331ad) records lmax16/24/32 expansion RMS.00359232/.000114168/9.61178e-6 with diagnostic coarse flags false. At lmax48, ntheta50, alpha.2 and up to3000iterations, the strict row passes: RMS9.92368e-8, area relative error1.75415e-14, sampled off-grid shape relative maximum3.34866e-7, runtime238.94s. The row and full four-order aggregate pass, and actual cycle0/MeshBlock-cycles0 is verified from the log. This removes the gross alpha1 oscillation without changing the physical seed or its acceptance thresholds. Fixed-l quadrature/tolerance refinement and actual solved-binary validation/enclosure remain separate requirements. The source-only known-ellipsoid expansion control F=(x-center)^T[I+(Gamma²-1)vhat vhat^T](x-center)-rh², with analytic gradient/Hessian, is recommended as an additional cheap diagnostic if future flow failures need to be separated from geometry; the reviewer has not run it.

For future axis-regular basis work, the existing map gives a=tanh(X/2)=(A+1)/2, t=a², eta=cos(R)=-2B/(1+B²), and q=a sin(R)=[(A+1)/2]*(1-B²)/(1+B²). Thus the proposed factor q^|m| has no square-root branch in computational A,B. Expanding the remaining w_m in t,eta enforces even-transverse parity, whereas merely multiplying an arbitrary Chebyshev(A,B) function by q^m does not. A natural new t/eta grid requires separate conditioning and near-hole-resolution assessment; retaining the old A nodes after squaring would produce highly clustered t nodes. Cartesian evaluation can represent q^m exp(im phi) as [(y+iz)/(b*(xi+1))]^m to avoid axis angular division, with xi=cosh(X). These are future design recommendations, not implemented repairs or new validation evidence.

A conservative represented-surface enclosure check can be made without another elliptic solve. For orthonormal real harmonic coefficients grouped by l, L_R=sum_l ||a_l||₂ sqrt[l(l+1)(2l+1)/(4pi)] bounds the spherical gradient magnitude globally, by Cauchy–Schwarz and the gradient addition theorem. If a dense angular sample set has covering radius delta, min R_true >= min R_samples-L_R*delta. For tensor-product polar/azimuthal nodes a conservative delta is max(first polar angle,pi-last polar angle,half the maximum polar gap)+pi/Nphi. Require |c_AH-c_hole|+inner_max strictly below that lower bound, with additional angular-refinement/shape-error margin. The same bound can certify positive radius between nodes. This certifies support containment in the represented finite harmonic shape; physical AH validity still requires the independent expansion and refinement checks, and does not resolve the solved-axis metric limitation. The parent has been given this optional numerical-certification method; it has not yet been implemented or run by the reviewer.

### Loaded-image provenance correction

The parent confirms the shared macOS install-name defect above and invalidates the interrupted fresh revalidation that actually loaded the old image. No acceptance is inferred from those results. The new supported workflow launches sequential fresh execed workers, each loading exactly one HiSpID image; the coordinator loads none. dladdr verifies the actual symbol image resolves to the requested file. Backend records a process-wide registry keyed by default_config symbol address, with the first loaded resolved path/disk SHA, rejects a different current hash for that address, freezes the Backend SHA and rejects later on-disk modification. This closes the identified same-address/rebuilt-path gap in the supported module/process workflow. Controlled fresh-process/path/disk checks are not a general claim that an arbitrary currently mapped binary can be hashed by reading a pathname. The separate-process proof remains pending at this writing. The exact single-build AthenaK controls remain valid and are not comparisons of two library images.

### Smallest proposed consistent C2 axis repair

The user/parent now requires resolving the axis limitation. Read-only inspection covers HiSpID_solver.cpp, HiSpID_spectral.hpp, TP_CoordTransf.c and the FD/ILU preconditioner. The current unknowns are4*N nodal V values, u=W+(A-1)V0, b=(A-1)Vvector, on Gauss-Chebyshev A/B nodes and periodic phi. A viable minimal change keeps this grid, mapping, nodal-value semantics, coupled physical collocation equations, cache layout and approximate FD preconditioner, but replaces each Fourier mode's meridional interpolant with an endpoint-constrained cardinal extension. It does not project away PDE residual components or introduce a rank-deficient unknown projection. A remapped t/eta regular basis is a stronger C-infinity alternative, but changes grid resolution and mode weighting; naive division by q^|m| of physical nodal Fourier coefficients can catastrophically amplify near-axis roundoff.

Let a=(A+1)/2 and sigma=sin(R)=(1-B²)/(1+B²). Each Cartesian correction component has the following sufficient **C2** mode conditions:

| Fourier m | At A=-1, on V | At each B=s=±1, on V |
|---|---|---|
|0|2 V_A-V=0|V_B=0|
|1|V=0; V_AA-V_A=0|V=0; V_BB+s V_B=0|
|2|V=0; V_A=0|V=0; V_B=0|
|>=3|V=0; V_A=0; V_AA=0|V=0; V_B=0; V_BB=0|

These follow by requiring u_m to have admissible constant/linear/quadratic Cartesian Taylor terms and O(rho³) remainder. The A Robin terms include u=(A-1)V: v_a=2V-2V_a and v_aa=4V_a-2V_aa at a0. At B=s, sigma_B=-s and sigma_BB=1; the m1 condition removes an inadmissible quadratic transverse term. Equivalently R_B=1 and R_BB=-s at the endpoint, so v_RR=v_BB+s v_B. All m>=3 may use O(rho³) with arbitrary smooth angular factor: this is C2, and not generally C3/C-infinity. It suffices for continuous physical H (metric second derivatives) and M (K first derivatives) if u,b are C2, with K from first derivatives of b. A classical smoothness claim beyond this requires stronger mode factors/conditions.

A and B extension operators commute inside each mode, so the corner conditions are compatible. Scalar corner values are not forced to zero; m>=1 leading factors vanish on both axis lines, with m1~a sigma~rho. Puncture corners themselves remain removed/excluded. This repair proves C2 away from punctures, not full Cartesian regularity of the scalar correction at the puncture; bounded puncture corrections may have a cusp. Corrected V remains finite at A=1, preserving u,b=O(1/r), the fixed W reference, and current asymptotic charge conditions. Physical collocation locations and puncture resolution are unchanged.

For an N-node Gauss-Chebyshev line, Q=T_N vanishes at every collocation node. Let P be its ordinary degree<N interpolant and choose lifting functions chi_j=Q*T_Kj. The corrected interpolant is P+sum_j c_j chi_j. For the appropriate endpoint functionals L above, C_rj=L_r[chi_j], c=-C^-1 L[P]. There are k=1/2/2/3 conditions on A and2k on B. The nodal values remain identical, so no new DOFs or PDE equations are introduced.

**Do not use only K=0..k-1 as the lift.** Independent tiny NumPy linear algebra finds severe growing conditioning even after scaling derivative rows: at N160 high-mode condition numbers are about5.86e9 for A and3.67e8 for B. Such a minimum-degree lift can generate huge between-node amplitudes. Frequency-spread lifts avoid this. For A use the first k of K=[0,N/2,N] (integer-adjusted); for B use paired even/odd K around0,N/2,N, i.e [0,1,N/2,N/2+1,N,N+1] for even compatible N. Use distinct parity pairs for odd N. Scale value/first/second functionals by1,N²,N⁴, respectively (m0 Robin/first byN²; m1 second Robin byN⁴). Tested N8/16/32/64/160/256 gives scaled conditions<=319, with high-mode N160 about271(A)/273(B), rather than growing powers of N. This raises the represented degree to about2N, while remaining rank3/6; it can be stored and evaluated as a low-rank lift rather than a complete enlarged tensor cube. Avoid K=[0,N,2N]: Q*(1+T_2N)=2Q³ is invisible to values and first/second derivatives at all collocation nodes, undermining collocation control. Independent off-grid validation is essential for any lift because its degree exceeds the number of collocation nodes.

Let F be nodal-to-Chebyshev coefficient conversion, L the scaled endpoint matrix, and H=-C^-1 L F. At the Gauss nodes chi=0 exactly,

D_new=D_raw+B1 H,
D2_new=D2_raw+B2 H,
B1_ij=T_N'(x_i) T_Kj(x_i),
B2_ij=T_N''(x_i) T_Kj(x_i)+2 T_N'(x_i) T_Kj'(x_i).

Build D2 directly: **D_new*D_new is wrong**, because differentiating the corrected function does not preserve its original homogeneous endpoint conditions, and a second application would reinterpolate into the wrong space. Mixed A/B derivatives are D_A^m D_B^m, since the tensor extensions commute; phi differentiation commutes with the same-frequency class. Keep the separate cosine-Nyquist Fourier D2. Only four mode classes0/1/2/>=3 are needed. The existing transform's (A-1) product derivatives and A/B-to-Cartesian chain rules remain valid when fed these actual corrected derivatives. Its FD/ILU stencil remains an approximation on identical nodal unknowns/grid; Krylov quality must be measured, not assumed.

Independent tiny matrix tests for N16/32/64/160 find the corrected derivative infinity norms remain bounded relative to the inherited matrices: the worst high-mode A D/D2 ratios are9.69/9.56 at N16 and5.63/5.57 at N160; B is9.23/7.74 at N16 and5.07/3.65 at N160. Scaled endpoint operator identities have defects<=8.5e-15. These are small dense linear-algebra controls, not native solver runs. Numerical endpoint evaluation must handle floating-point cancellation carefully: the mathematical zero/Robin identities can be imposed exactly for axis Taylor evaluation instead of trusting cancellation of large endpoint derivatives. Existing Chebyshev sampling halves its stored constant-mode basis; chi_K uses ordinary T_K with T0=1 and must use a consistent normalization.

The sampler must evaluate the SAME corrected tensor interpolant as the collocation operator; calling SpecCoef on unchanged nodal values alone reconstructs the original unconstrained interpolant. Apply/store both modal low-rank lifts, including mixed corner contributions. Axis values/gradients must use its unique Cartesian limits, not the old uncorrected zero-mode central value. Current outer/interaxis inverses and constant lab-frame pullback remain valid. An explicit new unknown-parameterization identifier and fresh binary SHA, solve/physical/charge/covariance gates are required; this is a mathematical solution change, not an API-only migration of accepted coefficients.

Useful manufactured controls are exact low-degree A polynomials f_m=(1+a)*a^m*(1+a²), m0..3, and B functions h0=1, h1=.5(1-B²)*(3-B²), h2=(1-B²)², hhigh=(1-B²)³, tensored with cosine/sine modes, including m4. They satisfy the stated C2 conditions exactly; the independent N16..160 controls require lift coefficients<=2.4e-14. A high-mode rho³ fixture is intentionally only C2, so fourth-order Cartesian FD convergence must not be required there. Smooth independent Cartesian fixtures can instead use V_m=(1+a)q^m*w(t,eta), with w=1+.1t+.07eta+.03t eta, leading to u=-2(1-t) Re/Im[(y+iz)/(b*(xi+1))]^m w, xi=coshX. Test m0,1,2,3,4,6 on both axis segments, generic points, translated/rotated frames and a Nyquist mode. Compare physical values/first/second derivatives with this independent regular expression and retain former failing solved-axis h sequences. Only root owns native builds/solves.

### Draft axis-header audit and independent precision control

The reviewer reads the drafted src/HiSpID_axis.hpp and the subsequent root-owned collocation/sampler integration without editing either. The half-T0 coefficient normalization is correct: the legacy stored Chebyshev constant coefficient multiplies T0/2, so its endpoint-functional column carries .5; a cardinal/nodal column does not. The orthonormal real Fourier matrix has amplitudes1/sqrt(Nphi) for m0 and the cosine Nyquist, sqrt(2/Nphi) for the other cosine/sine pairs. Its modal ordering matches the inherited real transform. The first Nyquist derivative vanishes on the Fourier nodes, while its separately constructed second derivative is -m² times the mode, and the off-grid coefficient sampler differentiates it normally. Mixed meridional/azimuthal derivatives commute because both sine and cosine of a given frequency use the same A/B extension. Direct modal D2, rather than squaring D, is required and is used.

At alpha_j=pi(j+.5)/N, x_j=-cos(alpha_j), Q=T_N:

Q'(x_j)=(-1)^(N+j+1) N/sin(alpha_j),
Q''(x_j)=x_j Q'(x_j)/sin²(alpha_j).

These signs and exact-node lift formulas in the header are correct. The latter second derivative follows from the Chebyshev ODE at Q=0. B lift frequencies occur in even/odd pairs, including odd N, and the selected endpoint constraints have the stated parity. The two tensor lifts automatically include corner terms in the sampler's product of corrected bases.

The proposed direct nodal setup is also correct. For each endpoint e=±1, the cardinal polynomial is ell_j(x)=Q(x)/[Q'(x_j)(x-x_j)]. Its endpoint jets are Q/(q_j delta), (Q'-Q/delta)/(q_j delta), and (Q''-2Q'/delta+2Q/delta²)/(q_j delta). Compute delta=-2sin²(alpha_j/2) at e=-1, or +2cos²(alpha_j/2) at e=+1. This avoids subtracting an almost-endpoint cosine and avoids multiplying an endpoint coefficient functional by a trigonometric forward matrix. It constructs exactly the same corrected interpolant in real arithmetic.

Independent small NumPy calculations reproduce the setup and compatible polynomial controls; they are not native builds or PDE solves. With plain double endpoint-functional setup, normalized maximum D2 error at N8/16/32/64/160 is approximately3.3e-12/6.5e-11/4.9e-10/1.65e-8/5.39e-6. The raw D2 on the same fixtures is8e-14/2e-12/3e-11/3.2e-10/3.6e-8. Rationalizing pairwise nodal differences as2sin[(alpha_i+alpha_j)/2]sin[(alpha_i-alpha_j)/2], and using raw barycentric off-diagonal D2_ij=2D_ij[D_ii-1/(x_i-x_j)] with row-sum diagonal, barely changes the lifted error. This identifies the low-rank setup/endpoint-functional cancellation as the dominant avoidable error; it does not justify relaxing the test solely on raw-D2 roundoff grounds.

A separate Python stdlib Decimal control uses55-digit arithmetic, a55-digit pi, Taylor sin/cos, endpoint Chebyshev integer jets, the analytic cardinal expressions, a tiny pivoted elimination, and direct barycentric D2. It tests N160 A classes0/3 and B classes0/3 against the compatible polynomials above. No external high-precision dependency or source mutation is used. Exact high-precision polynomial reproduction has normalized errors below3e-44. Rounding only the final matrices to double and performing ordinary NumPy matrix-vector products gives, respectively,2.75e-8,6.00e-9,3.96e-9,1.31e-12. Keeping high-precision matrices and rounding only nodal polynomial data gives5.54e-9,1.02e-9, effectively zero for B0's constant, and8.45e-13. Therefore the observed few-e-6 setup error is avoidable; N^4 sensitivity still limits the final double computation. The Decimal oracle also computes trigonometric setup quantities at high precision, so a higher-precision small solve alone does not guarantee all of its improvement if its RHS already contains double cancellation. On ARM/macOS verify numeric_limits<long double>::digits before assuming long double supplies extra mantissa bits. The parent plans compensated double-double setup/accumulation with double runtime storage; its code and native controls remain pending.

There is a cancellation-free endpoint-moment alternative, mathematically equivalent to the header's conditions. Write z_k=(K_k/N)² and d_k=e^(N+K_k)c_k at endpoint e; moments M_r=sum_k z_k^r d_k. For class>=3, the value/first/second conditions give

M0=-P(e),
M1=P(e)-e P'(e)/N²,
M2=-3P''(e)/N⁴-(1-1/N²)M0-(6-1/N²)M1.

Solve the small Vandermonde rows [1,z_k,z_k²]. Class2 uses only M0/M1. For class1, both A=-1 and B=e conditions are P_corrected''+e P_corrected'=0, so use rows [1,G(z_k)], where G(z)=z²+(6+2/N²)z, with RHS

M0=-P(e),
G_rhs=-3[P''(e)+e P'(e)]/N⁴-(1+2/N²)M0.

For a coefficient basis P=f_i T_i with f_0=.5 and f_i=1 otherwise, let p=f_i e^i, r=i²/N², and calculate1-r=(N²-i²)/N² without subtraction of nearby ratios. The high-class coefficient RHS is [-p, p(1-r), -p(1-r)(5-r)]; class2 uses its first two entries; class1 uses [-p,p(1-r)(1+r+2/N²)]. These formulas avoid cancellation of the common endpoint derivative pieces.

For a cardinal column, put tau=e*(e-x_j)>0, Z=N²tau, ell=Q(e)/[Q'(x_j)*e*tau]. The high-class nodal RHS is [-ell,ell/Z,ell*(1/(N²Z)-6/Z²)]. Class2 uses the first two entries. Class1 uses [-ell,ell*((6+3/N²)/Z-6/Z²)]. These are exact algebraic simplifications of the same endpoint jets. A moment solve maps d to c by the sign e^(N+K). At B, separate lifting columns according to total parity N+K: even/odd moment RHS is(M_plus±M_minus)/2; solve one3×3 or2×2 system per parity, with c=d at e=+1. No new unknowns, endpoint constraints or interpolation space are introduced.

For completeness, A class0 has only K0 and c=-p/e^N*(i²+.5)/(N²+.5) for a coefficient column, or c=-ell/e^N*(1-1/Z+.5/N²)/(1+.5/N²) for a cardinal column. B class0 solves the even/odd derivative rows with weight1+z_k and endpoint RHS -eP'/N², namely -p*r for a coefficient column or -ell*(1-1/Z) for a cardinal column. A constant B coefficient therefore receives exactly zero lift.

An independent double test of high-class endpoint-moment setup at N160 gives normalized compatible D2 errors1.41e-7(A) and1.76e-7(B), versus the few-e-6 derivative-functional setup. The moment Vandermonde's condition number is20.50 for K=[0,N/2,N], compared with about273 for the scaled derivative-functional system. These are setup comparisons, not proof of solved physical accuracy. Both the sampler and collocation constructor must use the same moment-equivalent extension if this alternative is adopted.

The root-owned integrated sampler's interaxis/outer-axis Cartesian first-gradient formulas agree with the independently derived limits above. Only m0 contributes the axis value/axial derivative, and m1 contributes transverse first derivatives. Sampling phi0 and phi=pi/2 correctly separates cosine/sine m1 after endpoint value constraints. Its rho<1e-10*b band uses an axis Taylor value through first order; that is an explicit finite band approximation, not an exact off-axis expression. Its omitted value is O(rho²) and omitted gradient O(rho) away from a focus, with geometry-dependent constants. Puncture foci remain excluded. Full contexts now retain six modal work arrays (192B/point for4fields), while sampler-only contexts retain compressed lift matrices. Keep those allocations in the conservative budget accounting. The existing large full-context allowance appears to have headroom, but its description should include the new arrays explicitly.

Fresh build/SHA, native tests, off-grid/axis h studies, binary convergence, charges/covariance and horizon checks remain required. A C2 high-mode rho³ remainder gives generally O(h) finite-difference truncation at an axis; requiring fourth-order convergence there would incorrectly reject the stated function class. Fourth-order checks should use a smooth Cartesian manufactured fixture. Absence of the former h^-2 discontinuity alone is necessary but not sufficient for an exterior-constraint acceptance claim.

The parent then implements AxisWide compensated setup arithmetic. Read-only source inspection confirms its endpoint integer jets, common row scales, fma product residual, three-step quotient refinement, two-component pi, sine reduction/Taylor20, direct-node Chebyshev recurrence and wide low-rank accumulation are consistent for the bounded nonzero initialization operands. The raw ordinary diagonal D_ii=x_i/[2sin²(alpha_i)] and direct off-diagonal D2 formula are correct. The build uses ordinary -O3, without fast-math reassociation; preserve that requirement for compensated arithmetic. The initialization sine arguments lie within[-pi,pi], and Taylor truncation after reduction to[-pi/2,pi/2] is well below the compensated mantissa error.

An independent Python emulation uses Decimal to emulate the exact fma residual and compare the same two-component operations with70-digit trigonometric arithmetic. At representative first/quarter/mid/last nodes for N4/7/160/256, including half angles and cosine arguments, the compensated pi error is2.995e-33 and the largest sine/cosine error is4.896e-32. This verifies the arithmetic route independently; it is not a compilation or execution of the native header. Root reports the unchanged standalone native polynomial gate now passes with worst normalized error2.33149e-7 at N160, endpoint defect8.13e-15 and nodal reconstruction1.00e-13. The reviewer has not independently run that executable. Runtime matrices, Fourier transforms, PDE arithmetic and sampler basis evaluation remain doubles.

A dedicated approach-to-axis manufactured control is still needed, for example rho/b spanning1e-2 through1e-10 at fixed points on both axis segments, with m0/1/2/high/Nyquist combinations. The sampler's basis() still evaluates subtracting polynomial/lift pieces in doubles and imposes endpoint identities only at the exact endpoint. Outside its tiny Taylor band an eps-sized value/Robin error may be amplified by1/rho in Cartesian gradients; the exact-axis limit formula alone cannot bound that effect. Compare against smooth independent Cartesian values/gradients, include the band threshold, and report the attainable floor. If necessary, a stable endpoint Taylor/factored evaluation can encode A=-1 using w=A+1:

class0: V=c0*(1+w/2)+w²R(w),
class1: V=c1*(w+w²/2)+w³R(w),
class2: V=w²R(w),
class>=3: V=w³R(w).

These express exactly the same endpoint conditions. In particular u=(w-2)V has no linear term for m0 or quadratic term for m1. Near B=e let w=B-e: class0 V=c0+w²R, class1 V=c1*(w-e*w²/2)+w³R, class2 V=w²R, and class>=3 V=w³R. This is a numerical evaluation suggestion, not a further basis change, and it must represent the same full corrected polynomial. A finite C2 Taylor band may instead be used with explicit value/gradient truncation bounds and near-focus scaling; no such additional repair has been implemented by the reviewer.

### Elliptic conditioning evidence and subsequent factored-grid proposal

Passing endpoint/manufactured identities does not prove coercivity or good conditioning of the full collocated elliptic operator. Frequency-spread lifts add degree~2N content which is invisible to nodal values; the condition of their small endpoint system is not the condition of the PDE. Parent reports tiny generic coupled solves stalling even with larger Krylov budgets or a dense Newton step. These failures are retained and do not establish a working repaired binary solver.

The initial raw weighted BL scalar Jacobian condition3.84e15 and tiny eigenvalue5.4e-14 were initially interpreted as spurious null modes. That inference is **withdrawn**: the physical compactification and sin(alpha)^6 sin(beta)^6 weights create row norms spanning about1e-12 to58.6. Independent reading of axis-matrix-old.json, axis-matrix-c2.json and axis-matrix-minimal.json confirms the corrected row-equilibrated comparison: original SHA2dbf450d… has condition502.56; spread-lift SHA9fe471642… has10034.94; minimal-lift SHA0c1773288… has380789.34. All are full rank under this comparison. The spread basis worsens conditioning about20times and the minimal basis about758times; the misleading absolute raw-eigenvalue interpretation is not evidence of a true kernel. The parent owns all Jacobian/solve executions, using separate fresh processes.

A minimum-degree cardinal high-class lift is exactly (A+1)^3 times the degree<N interpolant of V_j/(A_j+1)^3, and similarly (1-B²)^3 times the degree<N interpolant of V_j/(1-B_j²)^3. These are degreesN+2 andN+5 and avoid needless doubled degree, but they can still magnify endpoint nodal noise; the observed conditioning is worse. Neither endpoint conditioning nor minimum polynomial degree is an adequate solver acceptance criterion.

The subsequent proposal replaces the meridional grid with Gauss nodes of t=a² in[0,1] and eta=cosR in[-1,1], and uses V_m=(1+a)q^r P_m(t,eta), q=a sqrt(1-eta²). Then v_m=(A-1)V_m=-2(1-t)q^r P_m, with v=u-W or a Cartesian vector correction. **Do not use r=min(m,3) for every m:** for even m>=4 it permits only odd transverse powers rho³,rho⁵,… in analytic P(t,eta), while a smooth Cartesian m4 term begins rho⁴. The quotient would require a square-root singularity, losing spectral convergence. A parity-preserving C2 cap is r=m for m<=3, r=3 for larger odd m, and r=4 for larger even m. Every smooth mode remains representable because q^(m-r)=[t(1-eta²)]^((m-r)/2) is polynomial. Higher odd/even modes may also represent only C2/C3 remainders; do not claim C-infinity from the capped representation.

The new Gauss grid changes focus distances from O(b/N⁴) on inherited A/B to O(b/N²) on t/eta, so attenuation-shell/puncture resolution and all physical convergence claims need new evidence. Reusing the old A nodes after squaring would instead create severely clustered t nodes. Large-frequency q^-m division is undesirable. Even capped q^-4 division of Fourier-transformed physical nodal V is unsafe at joint axes: at N160, q_min is approximately4.8e-5 and q^4 approximately5.4e-18. The parent therefore chooses **modal P as the primary unknown**, rather than computing it by dividing transformed physical nodal V. This is a new grid and unknown parameterization; old checkpoint arrays are not directly compatible or accepted warm starts.

For a weighted-cardinal basis S*interp(P), S D S^-1 alone omits the weight derivative. The correct matrices on physical nodal V would be D1w=S D S^-1+diag(S'/S) and D2w=S D2 S^-1+2diag(S'/S)(S D S^-1)+diag(S''/S). For primary modal P, directly compute S P, S'P+S P', and S''P+2S'P'+S P''. With S_t=(1+sqrt(t))*t^(r/2), p=r/2:

S_t'/S_t=[p+(p+.5)sqrt(t)]/[t(1+sqrt(t))],
S_t''/S_t=[p(p-1)+(p²-.25)sqrt(t)]/[t²(1+sqrt(t))].

For S_eta=(1-eta²)^(r/2): S_eta'/S_eta=-r eta/(1-eta²) and S_eta''/S_eta=[-r+r(r-1)eta²]/(1-eta²)². Ordinary Chebyshev differentiation on the normalized coordinate2t-1 has physical t first/second factors2/4. To reuse the old A/B Cartesian map transform, V_A=aV_t, V_AA=a²V_tt+.5V_t; eta_B=-2(1-B²)/(1+B²)², eta_BB=4B(3-B²)/(1+B²)³; V_AB=a eta_B V_teta, with the same factors for mixed phi derivatives. The tensor derivatives commute. The sampler must make only the two meridional coefficient transforms on a modal P array, or first inverse-Fourier-transform virtual P before using legacy SpecCoef; transforming an array whose third slot already denotes mode index as if it were physical phi is incorrect.

The flat scalar operator provides an independent analytic test and proposed modal preconditioner. Let h²=4t/(1-t)², s²=1-eta², Dphys=b²(h²+s²), c=-2(1-t)q^r and actual Fourier frequency m. Then

Laplacian[c P exp(im phi)]=(c/Dphys)*B_rm[P],

B_rm[P]=t(1-t)² P_tt
 +(1-t)*[(r+1)(1-t)-2t] P_t
 +(1-eta²) P_etaeta-2(r+1)eta P_eta
 +[-(r+1)(r+1-t)+(r²-m²)(1/h²+1/(1-eta²))] P.

For full r=m the apparent axis singular terms cancel exactly; the capped representation retains the nonpositive potential proportional to r²-m². Dropping that term would solve a different operator. For internal GMRES RHS in weighted conformal FH/FM, Fourier transform first and solve B_rm P=[Dphys/(weight*c)] RHS_mode; c is negative. No physical-H factor belongs to internal FH. Only a RHS starting from physical H requires -psi^5/8, and physical M requires psi^10 plus the lab-to-local frame pullback. Replacing the flat vector Laplacian plus one-third grad(div) by independent scalar Laplacians is an explicitly approximate preconditioner; uniformly multiplying all vector components by4/3 also changes transverse components and is not the exact longitudinal/transverse operator.

Primary modal P removes the nodal V/S roundoff leak, but dividing Fourier-transformed physical residual modes by q^r in a preconditioner can still amplify low-mode numerical noise. At joint axes test a known low-mode manufactured residual for spurious high-mode preconditioned output. Any fixed invertible modal row scaling should be documented with residual/physical diagnostics; do not silently threshold or truncate modes to hide amplification.

The factored representation's exact axis values/first derivatives have a simpler independent form. For |x|<b, t=0 and eta=x/b: v0=-2P0, v_x=-2P0_eta/b, v_y=-P_cos1/b, v_z=-P_sin1/b. For |x|>b put D=|x|+b, t=(|x|-b)/D, eta=sign(x): v0=-4bP0/D, v_x=-2*(2sign(x)b/D²)*[(1-t)P0_t-P0], v_y=-4bP_cos1/D² and v_z=-4bP_sin1/D². All higher classes have zero first-axis derivatives; add W and its analytic gradient only for the scalar correction. These formulas avoid sigma/sinh division, and can validate the general coordinate-transformed sampler. The factored-grid implementation and all native/scientific gates remain pending at this writing.

The proposed flat B_rm operator has a useful continuum stability identity. It is formally self-adjoint with measure mu=t^r(1-eta²)^r dt deta: radial divergence is mu^-1*d_t[mu*t(1-t)²*d_tP], angular divergence is mu^-1*d_eta[mu*(1-eta²)*d_etaP]. Bounded P has zero natural endpoint flux. Integration by parts gives the negative of the integral of

mu*[t(1-t)² P_t²+(1-eta²)P_eta²
 +(r+1)(r+1-t)P²+(m²-r²)*(1/h²+1/(1-eta²))*P²].

For the selected r<=m this is strictly negative for nonzero bounded P. It supplies a sign/energy test, not proof that an arbitrary collocation matrix inherits coercivity. A conservative nonuniform finite-volume approximation for the flat preconditioner can use radial flux t^(r+1)(1-t)² and angular flux(1-eta²)^(r+1), zero boundary flux, cell-width denominators and the negative pointwise potential. Its dual-cell/mu weighted matrix has symmetric flux couplings and negative energy, avoiding an unjustified reflected-P endpoint closure. Such a preconditioner remains approximate for the curved coupled operator.

An independent smooth Cartesian manufactured fixture for the capped basis is P_m=q^(m-r)=[t(1-eta²)]^((m-r)/2). It is polynomial because the cap preserves parity. It yields v=-4b Re/Im[(y+iz)^m]/D^(m+1), where D=(r_plus+r_minus+2b)/2=b*(xi+1). Test m0/1/2/3/4/5/6 plus a Nyquist mode and translated/rotated frames. Distance derivatives D_i=.5(n_plus_i+n_minus_i), D_ij=.5[(delta_ij-n_plus_i*n_plus_j)/r_plus+(delta_ij-n_minus_i*n_minus_j)/r_minus] give an exact independent value/gradient/Hessian oracle away from either focus. For R=Re/Im[(y+iz)^m] and p=m+1:

v_i=-4b*[R_i D^-p-p R D^(-p-1)D_i],
v_ij=-4b*[R_ij D^-p-p D^(-p-1)*(R_i D_j+R_j D_i+R D_ij)
 +p(p+1)R D^(-p-2)D_i D_j].

These avoid chart inversion, axis divisions and Fourier derivatives. They also test that B_rm[q^(m-r)] equals q^(m-r)B_mm[1]. The reviewer derives these formulas but does not run a native test or solve.

An independent tiny NumPy control then evaluates the capped B formula against the analytic Cartesian Laplacian of this fixture at b1.3, points(.2,.3,.4),(2.2,.2,-.1),(.2,.002,.001),(-2.2,.002,.001), and m0through6. The normalized maximum difference is9.751e-14. The Cartesian Laplacian uses Laplacian(R)=0, Laplacian(D)=1/r_plus+1/r_minus, and the distance gradients above; it does not call the native library or map derivative implementation. This separately verifies the AF-factor zero-order term and the r²-m² potential sign. It is a formula check, not a discretized solve or grid-convergence result.

### First implemented modal-P source and sampler controls

The parent implements modal_P_C2prolate_v1 and reports a formerly failing tiny generic n8 solve now converges in2Newton/33Krylov iterations, .00745s at the default1e-10 tolerance, alongside16native controls. This is encouraging nonlinear evidence, not a physical accuracy/convergence acceptance. The reviewer reads the new header/solver and finds the A/B derivative chain, logarithmic weight derivatives, all nine derivative slots, modal ordering, two meridional coefficient transforms, half-T0 constants, orthonormal sampling normalization, direct Cartesian first-gradient chain and capped potential consistent. The preconditioner's scalar row scale weight*mu*c/Dphys is negative and correct; its scalar potential*Dphys/mu term follows by factoring the same row. Curved inverse-trace and scalar potential are azimuthal averages, and the4/3 vector Laplacian is an explicitly approximate block. The alpha/beta FD factors are correct for t_alpha=sin(alpha)/2 and eta_beta=sin(beta). Reflections in these cosine angle coordinates express the bounded polynomial extension; they are not physical-t/eta Neumann conditions. Discrete energy/conditioning and RHS-leakage controls remain necessary.

The reviewer independently reads validation/regular_modes.json and its Cartesian oracle script: SHA5ae1a9ba2076bccf1e577f2d996e14b5c4d179ae3a240025d709d51676c48532,16²×16,14scalar rows m0through6 plus cosine Nyquist8. Points cover x=0,1.2,±5.4 for b3 andrho/b1e-2through1e-10, including the exact axes. Maximum scaled metric-gradient error is4.914e-16 and exact-axis FD Hamiltonian error is<=2.091e-9 at h.002. All row flags and aggregate pass. The earlier chart-gradient failure is separately retained; the stable inverse a=sinh(X)/(cosh(X)+1), sin(R) root, and direct Cartesian gradient resolve it. These results verify the reported finite set and amplitude1e-4. Near-axis high-m perturbations can be much smaller than the seed metric's roundoff, so tiny physical dgamma error alone does not tightly measure their relative gradient; the generic collocation AD derivative control complements it. Additional off-grid rho~b points, a vector correction/K/axis-momentum manufactured check, and rotated/translated frames are recommended. The reviewer runs no native context or solver.

A private test-only translation unit can include the solver implementation to access its anonymous-namespace preconditioner, without adding a public production ABI hook. On exact BL (mu1,zero scalar potential), construct a pure orthonormal-m0 constant-P RHS as weight*c0/Dphys*B00[1], inverse-transform it into physical phi, and apply the modal preconditioner. Report spurious high-mode modal P and reconstructed off-grid field/gradient/Hessian errors; nodal q^r suppression alone can hide amplified coefficient tails. A uniform Fourier mode should not produce material high-mode physical ringing from row division. This check is proposed and pending.

For a vector oracle with u0=0 and only b^x=v on the same flat BL seed, compute Cartesian derivatives of the distance fixture above. Then Lxx=4v_x/3,Lxy=v_y,Lxz=v_z,Lyy=Lzz=-2v_x/3,Lyz0, and physical Kij=psiBL^-2 Lij. The independent expected constraints are H=-psiBL^-12*[(8/3)v_x²+2v_y²+2v_z²] and contravariant M=psiBL^-10*(Laplacian(v)+v_xx/3,v_xy/3,v_xz/3). Mode0 has nonzero axis momentum and supplies a useful Hessian check; mode1 tests off-axis K. These formulas require sampling/FD only, not a nonlinear solve. Current exporter/parser recognize the exact native parameterization; historical V arrays remain incompatible with new P arrays.

### First modal-P binary series fails physical acceptance

Fresh moderate_c2prolate records are read from validation/results.json, explicitly SHA5ae1a9ba…/modal_P_C2prolate_v1, grids24²×12,40²×20,56²×28 and the same generic unequal-mass moderate free data. All internal solves converge (254/414/553Krylov iterations,1.27/12.18/62.85s), with fine physical-equivalent g1 collocation extrema about1.6e-14(H) and2.8e-15(M components). Independent off-grid physical checks nevertheless **fail**: fine near H/M RMS3.2916e-4/.0037535 and bulk H/M RMS1.69196e-4/.0275579. Changing the FD step by factors2/.5 leaves the large bulk residual essentially unchanged, so it is not the calibrated small-step noise floor. Fine charges include Jy=-.02804,Jz=-.99966, compared with large changing coarse values and the earlier moderate result near.116/.338; charge stability fails. These are failed current-basis results, not acceptance of the repaired solver. The historical moderate_far0_regular.log name is not current modal-basis evidence; its timestamp/rows belong to the older series.

Read-only NumPy inspection of retained unknown arrays shows the largest m4 primary vector P at the joint-axis corner i=j=0 has magnitude approximately3354/4908/4759 at24/40/56, with alternating signs, and m2 maxima44.5/266.7/179.9. Such large P values are suppressed by q^r at the nodes but can create derivative/off-grid tails. They are not alone proof of floating RHS leakage: sharp puncture/attenuation structure on the new less-clustered grid can require large corner coefficients. The parent's private random-P collocation/sampler value+Cartesian-gradient comparison agrees to5.8e-23, reported as evidence against a reconstruction mismatch. A private low-mode RHS test reports high-P leakage near7e-6 at N64; its reconstructed physical/off-grid impact remains pending. ILU is a single approximate application, so the expected low-mode response is the actual ILU response, not an exact B solve.

The parent proposes analytic coordinate clustering while retaining the factored modal class: t=lambda*s/[1-(1-lambda)s], s=(1+z)/2 on a standard Gauss-Chebyshev z grid, and eta=tanh(k*zeta)/tanh(k) on a standard Gauss-Chebyshev zeta grid. Fixed lambda=.2,k2 give nonzero endpoint Jacobians .2(radial) and approximately.14657(angular), reducing the joint-focus nearest distance by about5.8 relative to the unclustered new grid. Both maps are monotone analytic on their closed real intervals and have analytic inverses there; C2 parity and smooth-mode representability are preserved. A new v2 parameterization fingerprint must identify the maps/constants, and all accepted results require fresh validation.

Let d=1-lambda,T=tanh(k). Inverses and derivatives are

s=t/(lambda+d*t), s_t=lambda/(lambda+d*t)², s_tt=-2lambda*d/(lambda+d*t)³;
zeta=atanh(T*eta)/k, zeta_eta=T/[k(1-T²eta²)], zeta_etaeta=2T³eta/[k(1-T²eta²)²].

The normalized radial polynomial coordinate z=2s-1 has z_A=2a*s_t and z_AA=s_t+2t*s_tt. Replace the unclustered header's hardcoded second-map1 in the P derivative term. The angular chain is zeta_B=zeta_eta*eta_B and zeta_BB=zeta_etaeta*eta_B²+zeta_eta*eta_BB. The sampler returns P_t=2s_t P_z and P_eta=zeta_eta P_zeta. Weight logarithmic derivatives in physical A/B remain unchanged.

For the FD preconditioner in uniform angles, t_s=lambda/[1-d*s]², t_ss=2lambda*d/[1-d*s]³;
t_alpha=t_s*sin(alpha)/2, t_alphaalpha=t_ss*sin²(alpha)/4+t_s*cos(alpha)/2;
eta_zeta=k*sech²(k*zeta)/T, eta_zetazeta=-2k*tanh(k*zeta)*eta_zeta;
eta_beta=eta_zeta*sin(beta), eta_betabeta=eta_zetazeta*sin²(beta)+eta_zeta*cos(beta).

For the same physical B coefficients ctt,ct,cee,ce, the angle-second/first coefficients are ctt/t_alpha², ct/t_alpha-ctt*t_alphaalpha/t_alpha³, and cee/eta_beta², ce/eta_beta-cee*eta_betabeta/eta_beta³. The capped B potential and RHS sign do not change. Cosine-grid reflection remains the coordinate-polynomial extension.

Manufactured exactness tests must match the new polynomial coordinates. P polynomial in physical t/eta is now rational/tanh in s/zeta and is not exactly represented at N12/13. Use P polynomial in actual s/zeta for exact chart-AD controls, independently differentiating the inverse maps in the oracle. The high-mode smooth Cartesian fixture P=q^(m-r) instead requires an interpolation-resolution study (or sufficiently high N), while constant-P modesm<=4 remain exact. Do not relax a derivative test to hide this distinction. The clustered proposal's implementation/physical checks remain pending at this writing.

The parent subsequently implements modal_P_C2prolate_mapped_v2. Read-only inspection of the header, sampler, preconditioner and revised native AD fixture finds z_A/z_AA, zeta_B/zeta_BB, physical eta/B weight logarithmic derivatives, FD chains and sampler inverse-map derivatives agree with the formulas. The AD fixture now uses a polynomial in actual s/zeta and differentiates the inverse maps independently. No mapping sign or missing-factor defect is found. The fixed constants .2/2 and new native parameterization name are explicit. Fresh build/provenance, numerical controls and physical binary acceptance remain pending; this source audit does not promote the failed v1 series.

One metadata refinement is reported: current JSON stores the descriptive parameterization string, which initially remains identical for v1/v2, rather than the exact native ID. The SHA check protects ordinary replay, but the description alone cannot fingerprint the changed grid in typed migrations or human review. Store/enforce a separate exact unknown_parameterization_id and explicit mapped constants, and retain accurate descriptions for archived modal_v1 images. Existing v1 records must remain historical and must not be relabeled as v2.

A mathematically equivalent preconditioner-only Fourier calculation can reduce avoidable mean leakage: use a compensated row mean, keep the normal normalized m0 sum, and project rhs_phi-mean onto each nonzero mode. Exact nonzero Fourier rows have zero sum, so this does not threshold/truncate or alter the real-arithmetic transform. Constant physical-phi rows then produce exactly zero nonzero modes before q^-r division. It can reduce dominant near-axis m0 roundoff; other low-mode numerical leakage and true aliasing still require the private/off-grid tests. This improvement is proposed, not independently executed by the reviewer.

The reviewer reads moderate_c2mapped results through80²×28, SHA d2de7b9fe0de…, with all acceptance flags still false. At N80 near H/M RMS2.1317e-5/.00024374 and bulk H/M RMS2.8667e-5/.00245264 remain outside the requested physical gate. Bulk momentum is nonmonotone over the last refinements. Jz at radii100/200/400 is.35194/.39787/.22647, rather than a smooth stable radial fit, so its extrapolation remains unaccepted. These radius oscillations are not explained solely by a simple linear-in-radius odd1/r tail; unresolved radial/corner structure remains a possible cause. Coordinate and manufactured-source audits do not establish a physical binary solution, and the failed v1/v2 numerical sequences must remain preserved.

The N80 mapped record is also **internally unfinished**, not a converged discretized solution: status1/Krylov iteration limit,3Newton/613Krylov steps, weighted M extrema1.93e-6/1.17e-6/1.20e-6 and unweighted M extrema.789/2.936/.774. Its physical errors/charges cannot be attributed solely to truncation of a converged spectral solution. The N24/40/56 mapped records do converge internally, but all fail independent physical acceptance. A fresh converged N80 result is required before treating that refinement as an accuracy trend.

Current replay metadata now includes exact unknown_parameterization_id and collocation_maps; checkpoints.py enforces both against the backend, and run_validation.py uses exact ID in its resume guard. This resolves the reported descriptive-string ambiguity without relabeling historical v1 arrays. The preconditioner implements compensated mean subtraction for nonzero Fourier rows, with m0=sum/sqrt(N), consistent with the equivalent transform proposed above.

Read-only inspection of mapped_cache_equivalence.json records serial derivative-cache locality and centered preconditioner projection between old mapped SHA d2de7b9fe0de67403a88ff520fa164ab4b3076c84d321c7c36180b892991f652 and new SHA722df13ac2a86cc390b0c7a677d2c2fea92034dd48d5e44f5b48100f08dc9fdd. Separate fresh processes/loaded-image checks are explicit; retained N80 off-axis gamma/K/psi/conformal data/corrections compare exactly, and residual/JVP absolute differences are zero. This supports the documented implementation-equivalence check, not numerical acceptance of the unfinished payload. Its five exact-axis points show H RMS about5.18e-6 and M RMS about7.50e-5 over h=.008,.004,.002,.001, with no former h^-2 divergence; the nonzero truncation/unfinished-solve errors still fail the physical gate. The reviewer neither loads a native context nor runs a solver for this inspection.

### Independent compactification and leading-angular-tail audit

The r3/4 cap permits high-m C2/C3 remainders, but does not intrinsically prescribe a wrong vacuum solution. Away from punctures the actual free data and elliptic coefficients are smooth; an exact positive-psi C2 solution of the uniformly elliptic coupled constraints would gain regularity. Finite collocation must enforce that property numerically. More decisively, the parent's retained v1 cutoff-mode test already reproduces bulk M RMS.027588 at cutoff4 versus.027558 for all modes. Modesm<=4 use fullr=m and are smooth on each non-focus axis segment. Thus capped modesm>4 are not needed for the dominant reported v1 momentum failure. This does not rule out finite meridional underresolution or inaccurate endpoint angular structure in low modes.

For far0, f,g approach1 at infinity. Individual seed vacuum identities make the leading r^-3 scalar and vector free-source terms cancel under linear superposition; trace-projection/metric-difference corrections enter at higher order. Let Q_u(n)=lim(R*u) and Q_b(n)=lim(R*b), with R the physical spherical radius, and use the mapped modal interpolant endpointt1:

Q(n)=-4b*sum_m[(1-eta²)^(r/2)*P_m(1,eta)*F_m(phi)], eta=n_x.

The continuum leading scalar equation is Delta_S2 Q_u=0, so bounded smooth Q_u is constant. For each Cartesian vector component, define the tangential gradient G_i=(delta_ij-n_i*n_j)*partial_nj and D=-n dot Q_b+G_j Q_b^j. The leading vector equation is

Delta_S2 Q_b+(1/3)*[-2n*D+G D]=0.

Its momentum/Kelvin kernel is Q_b=C*[7P+n(n dot P)], with arbitrary constant vector P and scalar normalization C. It has even parity under n→-n. Enforcing only O(1/R) decay does not directly enforce either angular equation in a finite polynomial approximation. Very small compactified residual rows can leave large angular defects. These are necessary asymptotic equations/diagnostics derived from the PDE, not additional published finite-radius boundary prescriptions. Do not impose b=O(R^-2) for generic net momentum: paper v3 §III.A PDF pp8–9 states only u,b→0 as physical boundary conditions; Fig1's quadratic b falloff is specifically associated with no net linear momentum. The same pages explicitly acknowledge logarithmic terms in unfiltered boosted superpositions and potentially algebraic spectral convergence, with attenuation-shell small scales limiting practical accuracy.

The reviewer computes a small independent NPZ/NumPy diagnostic without loading native code. Retained arrays are reshaped(mode,j,i,field). The exact radial endpoint cardinal weights for ascending Gauss nodesz_i=-cos(theta_i), theta_i=pi*(i+.5)/N, are

ell_i(1)=(-1)^(N-1-i)*sin(theta_i)/[N*(1-z_i)].

Collapse the radial nodal array with these weights, transform its remaining angular samples using the explicit Chebyshev cosine matrix, halve its zeroth coefficient, and evaluate at zeta=atanh[tanh(2)*eta]/2. Differentiate the polynomial and analytic inverse map directly. The Fourier factors retain native orthonormal real-mode normalization; all scalar/vector field factors are then multiplied by -4b=-12. This uses only saved coefficients and mathematical formulas, with single-thread NumPy runtime; no source/build/solver/sampler mutation or native execution occurs.

At mapped N56 (internally converged but physically failed), scalar Q anisotropy RMS is9.1994e-6 and vector odd-parity Q RMS1.01615e-3. The correction's leading H*R³ RMS is.06798079 and M*R³ RMS2.76583;160→256 Gauss-Legendre polar refinement changes these by7.8e-8 and5.2e-5 respectively. For the unfinished mapped N80, these leading norms are.0404695 and4.81169, with vector odd Q RMS approximately4.46e-4. The N80 angular norms still change modestly between160/256 polar points, so quote their shown accuracy rather than a tight gate. The large asymptotic angular derivatives coexist with small angular field amplitudes and explain why a tiny absolute far-field residual alone is weak evidence.

For reproducibility, in eta/phi coordinates write E=(1-eta²,-eta*n_y,-eta*n_z) and F=(0,-n_z/(1-eta²),n_y/(1-eta²)). A Cartesian component of beta=Q/R has R²*partial_j beta^i=-n_j Q^i+E_j Q_eta^i+F_j Q_phi^i. Apply the flat longitudinal operator to this tensor for the leading ADM integrand; alternatively Delta_S2 Q=(1-eta²)Q_etaeta-2eta Q_eta+Q_phiphi/(1-eta²), and use the D equation above for M*R³. Analytic Kelvin Q=7P+n(n dot P) gives maximum angular-M defect4.69e-15. A toroidal controlQ=S cross n gives J/R=-2S/3 within2.78e-16. These independent controls check signs, tensor indices and quadrature normalization.

The represented endpoint fields have tiny converged toroidal l1 content: direct leading J/R for N56 is(1.69e-10,-1.71663e-8,-6.52111e-8), and for unfinished N80 is(3.48e-11,-2.05706e-9,-5.28745e-9). The continuum integration-by-parts identity is J/R=-average_S2(n cross Q_b), verified against the direct derivative integrand at refined quadrature. At radii100–400 these converged coefficients are too small to explain the large finite-radius J oscillations by a true simple linear-in-R angular-momentum divergence. Nonphysical higher angular modes, finite-radius radial structure and quadrature aliasing remain distinct possibilities.

Crucially, fixed16/24/32 polar quadrature is not adequate for these leading angular tails. With40 azimuthal points, the direct N56 leading Jz/R sequence at polar16/24/32/48/64/96/160/256 is1.45870e-4,7.89633e-5,1.31826e-5,1.28439e-5,9.15283e-7,-6.99920e-8,-6.52111e-8,-6.52111e-8. For N80 it is1.52338e-5,-2.47251e-5,-1.40060e-5,7.68858e-7,2.32545e-6,7.89376e-8,-5.28749e-9,-5.28745e-9. Thus native/independent ADM agreement at identical16-point polar nodes can reproduce the same substantial angular aliasing. A finite-radius quadrature refinement at fixed radii and retained checkpoint is needed before attributing all measured J oscillation to the physical interpolant. This endpoint-only diagnostic does not compute full native finite-radius ADM charges or repair the already failed off-grid physical constraints.

### Exact modal FD block inverse and expanded Cartesian oracle

The parent replaces approximate ILU0 with direct block elimination of the same five-point modal FD preconditioner, preserving the actual spectral PDE, unknown parameterization and sampler. Read-only source audit of ModalBlock/Sparse finds the ordering/contractions correct. Write S_j=A_j-diag(lower_j)*T_(j-1), T_j=S_j^-1*diag(upper_j). The implementation forms each transfer column by a strided GSL vector view, subtracts the lower diagonal on the correct left row, saves independent pivot permutations, solves f_j=S_j^-1*(rhs_j-diag(lower_j)*f_(j-1)), and substitutes x_j=f_j-T_j*x_(j+1). Sharing sine/cosine partners and three vector components is valid for the azimuth-averaged normalized FD blocks; scalar potential retains a separate factor group. Each physical residual still gets its own negative weight*mu*c/Dphys row scale, and the4/3 vector approximation remains explicit.

The allocation screen explicitly adds2*(nphi/2+1)*nb*(16na²+24na) bytes on64-bit: two dense double matrices (LU/transfer), two diagonal coupling arrays, and size_t permutations for each scalar/vector/frequency group. This exactly counts those persistent arrays. Temporary factor RHS and one solve workspace are small relative to the existing conservative per-point stencil/Krylov/work-array allowance. Budget rejection occurs before allocation; the new dense factors must not be omitted from documented memory estimates. The reviewer does not build or execute this native implementation.

Parent-owned private tests assemble actual FD matrix-times-random-known-vector and invert it using the saved block factors, for N8/16/32/64, all four fields and real-mode partners. The reported maximum inverse errors are1.17e-15,6.05e-15,1.44e-14,2.45e-14. Centered purem0 residual tests give exactly zero spurious higher-P modes and zero reconstructed off-grid value/gradient leakage. Read-only inspection of the fixture notes that g=0,R=K=0 makes scalar/vector normalized blocks identical; an additional bounded variable scalar potential would meaningfully test their distinct grouping. This refinement is recommended, not falsely counted as already run.

The separate-process modal_block_equivalence.json binds old SHA722df13ac2a86cc390b0c7a677d2c2fea92034dd48d5e44f5b48100f08dc9fdd to new SHA25ca63784d1d18175bc85b06d20f940577bb6106d1d11792307cfaba1f63ce37. Controlled fresh-process/loaded-image checks are explicit. Retained mapped N80 fields and equation residual/JVP absolute differences are zero, establishing the stated operator/sampler equivalence for that payload. The17Python native controls also pass in the parent-owned log. These are implementation controls, not binary physical acceptance.

Read-only inspection of the expanded regular_modes_mapped.json reports28scalar/vector rows at64²×16, amplitude1e-4, mapped-v2/cache SHA722df13… . Cartesian distance-based values/gradients/Hessians remain independent of the map; points now include transverse radii. Allflags pass: largest correction error1.17e-18, scaled metric-gradient error5.62e-16, K tensor error1.39e-17, physical H oracle error3.69e-9 and M oracle error3.75e-16 at h=.002. Source audit confirms the bx oracle's Lxx=4v_x/3,Lxy=v_y,Lxz=v_z,Lyy=Lzz=-2v_x/3 and physical H/M formulas above. The numerical controls support the manufactured scalar/vector basis and sampling operators; they do not establish accuracy of unresolved binary free-source structure.

The parent subsequently obtains an internally converged cache N80 solve (5Newton/896Krylov,296.95s), then the block-preconditioned N80 solve (5Newton/64Krylov,31.33s). The latter log still gives near H/M RMS2.13270e-5/.000248068 and bulk H/M RMS2.86931e-5/.00243548, failing physical acceptance despite weighted extrema~1e-16. Its N104 refinement converges5Newton/54Krylov in42.43s, near H/M RMS5.71749e-6/9.49288e-5 and bulk3.97901e-6/.000231209. This improves the previously unresolved error but bulk M remains above the preliminary1e-4 gate and far above strict1e-6. The earlier unfinished N80 payload/endpoint diagnostic remains identified separately; it must not be relabeled as this converged solve. Current finer-resolution, quadrature, covariance and horizon evidence is still pending.

The parent then reports a N128 far40 control with bulk M RMS2.761e-4 versus far0's2.765e-4. This weakens the hypothesis that unfiltered asymptotic logarithms are the dominant cause of the current momentum failure; both remain failed constructions. The endpoint angular/quadrature defects above are diagnostics, not a claimed single explanation for all off-grid errors.

An important untested distinction is identified: the native chart AD fixture checks all nine A/B/phi derivative slots, and private random-P consistency checks Cartesian values/first gradients, but neither independently checks the inherited chart-to-Cartesian **Hessian** at joint-axis/focus nodes. The continuous distance oracle plus off-grid physical FD tests validates the sampler's first derivatives, not the collocated Cartesian second derivative used by the PDE. Inverting low Fourier modes into physical-phi rows before the Hessian transform can expose subtractive cancellation; subsequent division by q^4 in a modal preconditioner can magnify small phi-dependent errors into large corner P. This is a plausible mechanism pending direct tests, not a diagnosed defect.

The proposed independent control compares every transformed Cartesian Hessian slot with the distance-expression Hessian given above, including first/last meridional and ordinary nodes, m0through8 and scalar/bx components. A test-only flat cache with g0 or a flat Jacobian at zero corrections isolates FH=Delta u and FM=Delta b+(1/3)grad(div b), removing huge near-focus BL seed-LapPsi cancellation and nonlinear A² from that operator check. Report absolute and amplitude-normalized errors, with m0 first. Constant P for m<=4 is exactly represented on the mapped grid; m>4's smooth fixtureP=q^(m-r) needs sufficiently highN/interpolation refinement because the physical t/eta dependence is nonpolynomial in mapped s/zeta. Manufactured independent RHS→flat discrete JVP inversion should measure reconstructed Cartesian fields/Hessians as well as nodal P, and should not interpret a small internal residual as independent physical acceptance.

For purem0 the scalar oracle trace can be evaluated from analytic distance identities with phi-independent distances, instead of summing rounded Cartesian Hessian diagonals which may inject their own artificial Fourier leakage. Component Hessians remain useful independent comparisons. Near a focus Cartesian x-center cancellation creates a coordinate floor O(eps*b/r); either bound/report it or evaluate analytic distances r_plus/minus=b*(xi∓eta) and distance derivatives in extended precision. Ordinary points plus joint corners distinguish that input-coordinate floor from a Hessian-transform error. The reviewer sends these test-design cautions but executes no native context/job while the parent's fine solve is active.

If that control confirms a transform-cancellation defect, the same interpolant admits a direct Cartesian derivative evaluation, without another basis/PDE change:

v_m=-4b*G_(r,m)(y,z)*D^(-r-1)*P_m(t,eta),
D=(r_plus+r_minus+2b)/2, t=1-2b/D, eta=(r_minus-r_plus)/(2b).

For fullr=m<=4, G is exactly Re/Im[(y+iz)^m], giving polynomial Cartesian first/second derivatives without trigonometric cancellations. For higher capped modes G=rho^r*cos/sin(mphi) has first/second derivatives with explicit bounded rho^(r-1),rho^(r-2) factors away from the excluded focus and the exact axis limits above. Use distance derivatives D_i,D_ij already given, t_i=2bD_i/D², t_ij=2bD_ij/D²-4bD_iD_j/D³, and eta derivatives from the difference of the two distance gradients/Hessians. Chain discrete P_t/P_eta/P_tt/P_teta/P_etaeta through these coordinates and apply the product rule, then sum modes. Keeping low modes separate until after Cartesian differentiation may protect their regularity from physical-phi cancellation. This is a proposed stable evaluation route pending evidence; neither its native implementation nor a repaired numerical result is counted here.

The parent executes the independent Cartesian **operator** control on SHA25ca6378… . Read-only inspection of regular_operators_mapped.json and its source confirms28scalar/bx rows per grid, allcollocationpoints, m0/1/2/3/4/5/6/8 plus sine partners. A fixed zero-correction BL seed LapPsi floor is separately measured and subtracted only from the conformal mapped-operator discrepancy; full physical differences and the seed floor are preserved. N16 highm5/6/8 rows fail as expected from mapped-coordinate interpolation; N32 andN64 allpass, with largest normalized conformal errors5.37233e-9 and4.08085e-11, and full physical errors2.82614e-12 and1.46486e-14. The N64 seed conformal cancellation floor is2.16715e-5 near a focus, while the large puncture conformal factor suppresses its physical effect. The earlier seed-cancellation-confounded artifact is retained. These results weaken a gross transformed-operator defect for these manufactured fields; they do not exclude rounding amplified during a general solve.

Scope is stated precisely: scalar andbx controls independently constrain Lap(v),v_xx,v_xy,v_xz, but not individualv_yy,v_zz,v_yz. Additionalby/bz perturbations or direct private-TU comparisons of allsix transformed Hessian slots complete that coverage. For b^i=v, M=psi^-10*[ddv_i/3+e_i*trace(ddv)] and H=-psi^-12*[2|grad(v)|²+(2/3)v_i²]. This is a proposed scope refinement, not counted as already executed.

The parent next drafts meridional derivative rows in difference form sum_(j!=i)D_ij*(P_j-P_i), with the D2 diagonal replaced by minus its offdiagonal sum. Read-only source review confirms this annihilates constantP exactly and expresses the same analytic polynomial derivative. differentiate() is a distinct function from along(); AxisDerivatives uses it only for raw meridional D/D2/mixed operations, while coefficient and Fourier transforms keep their full matrix multiplication. Aliased input/output is not used in these calls. The native library remains25ca at this writing, so new rounding behavior has not yet been compiled or accepted. Any new derivative/operator evidence must bind the fresh SHA.

The private analytic flat RHS→GMRES audit uses constantPm0/m4, the independently derived B_rm source, and Cartesian distance-expression value/gradient error, with requested relative linear tolerance1e-13 and up to2000Krylov iterations. That strict target is useful as a conditioning/rounding diagnostic, but a global weighted residual is not a local joint-focus/high-mode accuracy norm. Record the analytic exact-P versus native JVP residual relative to RHS, the final explicitly reevaluated residual, and the reconstructed physical field error separately. A tolerance near the rounded operator's attainable floor may require many iterations without improving the physical representation. Preserve the old strict outcome and try the unchanged target after improved derivative summation before assigning a measured numerical floor; do not silently relax any binary physical acceptance criterion.

For the particular q^-4 leakage hypothesis, additionally probe modes1/2. Mean-centering annihilates purem0 physical-phi residuals exactly, and a purem4 RHS already carries q^4, so leakage into another cappedr4 mode does not acquire a large extra q-power ratio. A lower mode1/2 RHS spuriously projected onto r4 can acquire q^-3/q^-2 before inversion. Thus0/4 alone does not tightly bound the proposed low-mode-to-high-mode amplification. This extra control is recommended after the current strict audit, not falsely listed as executed.

The parent builds the row-difference derivative method as SHA c9186bc59fcd7abcd93989b14e380bdb69d40bbff2a2af15a1aec65de3fa8682. Read-only inspection of difference_derivative_equivalence.json confirms separate-process image verification, bitwise identical off-axis sampled fields, residual difference5.81e-16 and JVP difference2.43e-17. The mathematical interpolant is unchanged, but the rounded operators are explicitly **not** bitwise identical. At five solved N152 axis points on all three map-axis segments, momentum RMS approaches1.463e-5 as the Cartesian FD step decreases .008→.001; Hamiltonian RMS is about2.40e-7. The earlier h^-2 defect is absent, but this momentum value still fails the strict physical gate. The parent also reports a fresh default-map N80 solve retaining near/bulk momentum RMS2.48e-4/.002435. Improved differentiation summation therefore does not remove the binary's large physical error.

The regular-basis method now separates the primary modal unknowns P from physical-phi residual rows. The represented scalar correction and each Cartesian vector component are -2(1-t)q^r P times the orthonormal real Fourier basis, with q=sqrt(t)*sqrt(1-eta²). Frequencies0through4 use r=m; higher odd/even frequencies use3/4 respectively. The latter are C2 at ordinary axes rather than universally smooth, while every smooth Cartesian harmonic remains representable. Derivative construction multiplies P by the analytic weights before Fourier inversion; it never obtains high-mode P by dividing rounded physical-phi V by q^r. The coefficient sampler and collocation operator share the same two mapped Chebyshev interpolants and Fourier normalization. A correction with a capped high mode need not exhibit fourth-order FD convergence exactly on an axis; compatible smooth manufactured fixtures do. All accepted solved results still require axis approach/step and solver-resolution studies, not only these fixtures.

The exact modal block preconditioner is an inverse of the explicitly assembled five-point FD approximation to B_rm, not an exact inverse of the pseudospectral PDE. Factoring Delta(c P exp(i m phi))=(c/Dphys)B_rm P fixes its row scale: scalar weight*mu*c/Dphys and approximate-vector weight*(4/3)*mu*c/Dphys, where c=-2(1-t)q^r and mu is the azimuthal mean trace(h^-1)/3. The scalar derivative of the nonlinear algebraic source contributes potential*Dphys/mu to B_rm; no extra c is needed. Both radial and angular mapped first/second derivative chains are included. Real sine/cosine partners share one factor, the three vector components share another, and the scalar factor remains separate when its potential is nonzero. The dense block Schur update, strided transfer solves, forward/back substitution and mode/component indexing have been read-only audited. The known-vector matrix-times-inverse checks test that assembled FD inverse, while the independent analytic RHS-to-spectral-JVP tests below test its use within GMRES. These are distinct claims; neither substitutes for physical binary constraints.

The parent next builds a focused-map diagnostic with lambda=.05,kappa=3, identified through the new exact native string modal_P_C2prolate_map_v3_r<lambda>_k<kappa> and the collocation-map ABI query. The underlying formulas and endpoint regularity are unchanged: t_sigma(0)=lambda, t_sigma(1)=1/lambda, eta_zeta(±1)=kappa*sech²(kappa)/tanh(kappa). For equal large meridional N, the nearest joint-focus distance is approximately b*pi²*(lambda+eta_zeta(1))/(8N²). Thus .05/3 improves it by about4.35 over .2/2, but weakens the radial infinity spacing by a factor4. The map and its inverse are analytic and monotone on the real intervals. Their complex singularities move closer to the interval (the rational map pole lies at z=(1+lambda)/(1-lambda), and tanh's nearest pole at zeta=i*pi/(2kappa)); this can make mapped nonpolynomial fixtures converge more slowly. Pole locations describe an available analyticity region, not an unconditional error bound. Prolongation in raw Chebyshev coordinates is valid only within the same map; a cross-map migration would require explicit physical-coordinate interpolation and separate provenance. The native identifier, parser/checkpoint guards, Python map query and prolongation dispatch are read-only inspected. Historical v1 dispatch is also corrected.

Inspection of focus05-k3-build-tests.log confirms all nine weighted chart-jet values/derivatives agree with an independent map-aware AD fixture to2.053e-13, and the private random-field collocation/sampler Cartesian values/first gradients agree to5.82e-21. The initial hardcoded-default-map oracle failure is preserved separately. The strengthened flat inverse test now constructs every cached scalar/vector/L coefficient from exact constant delta jets, psi1 and zero free sources, preserving only the chosen row weights. This removes the previous BL quotient/connection rounding from the exact-flat control; the BL-rounded control is retained with its separate label. Analytic modes0/1/2/4 use the independent B_rm constant-P source and independently computed Cartesian distance-expression value/gradient at(.7,.8,.4). Exact-P native-action floor, final explicitly re-evaluated linear residual, physical value/gradient error and other-mode P amplitudes are printed separately.

The focused-map exact-flat and BL-rounded inverse rows through N32 all pass the declared relative1e-13 target. At N32, exact-flat m0 takes105Krylov iterations versus13 for the BL-rounded control; their sampled field/gradient errors are3.79e-17 and1.54e-16 respectively. These different iteration counts are not a continuum operator difference. At N64 the BL-rounded m0 row **fails** after2000iterations: exact-P linear floor6.85e-16, final relative residual3.899e-13, sampled value/gradient error1.837e-16 and other-mode P maximum1.690e-7. Its m1/m2/m4 rows pass, with m1 final relative residual6.09e-14 and sampled field/gradient error3.14e-17. The completed exact-flat N64 m0 row likewise **fails** after2000iterations: exact-P floor6.884e-16, final relative residual5.980e-13, sampled field/gradient error2.179e-16, other-mode P maximum3.897e-7. Exact-flat m1/m2/m4 pass; m1 takes303iterations, with final residual8.326e-14 and field/gradient error1.725e-16. The aggregate strict manufactured flag remains failed. Since the weighted RHS magnitude is small, this diagnostic does not imply that the production absolute1e-14 weighted target must fail, nor does small reconstructed error authorize weakening the physical acceptance gates. The exact-P floors near1e-15 are materially below the missed target; the larger attained residuals therefore reflect the iterative/rounded coupling, rather than disagreement of the analytic constant-P source with the discrete operator at that scale.

The m3/4 core ringing is not explained by the high-mode exponent cap: these frequencies retain full r=m. Large corner P values also do not by themselves prove roundoff leakage. Near a puncture D→2b, so if a physical component has Taylor coefficient v_m≈c_m*rho^m*cos/sin(mphi), its auxiliary coefficient approaches P_m=-2^(m-1)*b^m*c_m, up to the explicit Fourier normalization. For b3,m4, even a physical mode amplitude1e-4 at rho.09 corresponds to |P| about988; amplitude5e-4 corresponds to about4940, comparable with observed corner amplitudes. Narrow g/f shells can require rapidly varying but legitimate Taylor coefficients. With g0 and inner_flatten1, the corrections satisfy homogeneous flat Laplace/Navier equations in the core; bounded solutions extend analytically through the excluded puncture. Such Taylor solutions are representable by this basis. A focus-resolution experiment is consequently principled; it does not establish that underresolution is already the diagnosed cause.

If stronger mapping and anisotropic meridional resolution do not produce physical convergence, a principled fallback is a local puncture-core/exterior domain decomposition, keeping regular Cartesian modal factors and matching correction values and normal derivatives across smooth interfaces. It gives the narrow attenuation shells independent resolution without moving all exterior/infinity nodes inward. It is a larger implementation requiring interface-constraint and charge checks, not a recommendation to silently alter the current problem. Changing g widths/operator choice changes the PDE and needs its own declared configuration plus binary-horizon enclosure evidence. Removing capped-mode potentials, filtering solved m3/4 modes, or restoring the old axis-nonregular unknowns would not be an accuracy repair of this formulation.

The current acceptance gaps remain substantial. No repaired-basis generic moderate binary has passed the independent physical C/D gate; the N128 far0/far40 similarity and the fresh N80 row-difference result remain failed evidence. Fresh rotation/translation covariance, independent charge angular refinement before radial fits, exterior/attenuation-region separation, and solved-axis resolution checks are required for each accepted library/map. Earlier high-spin/high-boost numerical evidence belongs to historical libraries and bases; current chi.95/v.885 and combined binary range claims need fresh qualified solves. Exact-seed AthenaK/FastFlow controls validate geometry import and initial-time finder behavior, not binary horizons or enclosure of modified binary regions. Until actual solved binary surfaces with accepted expansion and refined shapes enclose the constraint-modification regions (g<1 and optional interior operator), an exterior-vacuum binary construction has not been demonstrated. The f/F free-data distinction is made explicitly below.

Two further discriminating diagnostics are proposed, without reviewer numerical jobs. First, the inherited row weight [sin(alpha)sin(beta)]^6 scales as N^-12 at a joint-focus Gauss corner. Its relation to physical focus distance changes when replacing the old A/B grid's O(N^-4) nearest distance with the regular mapped t/eta grid's O(N^-2) distance. A tiny absolute weighted residual can therefore hide large raw core rows; solver stopping/Krylov norms may neglect them while a global polynomial exports their interpolation error. In g0 with inner_flatten1 the correction operator is regular flat Laplace/Navier, so there is no nonzero singular free source requiring this strong suppression there. Record raw residuals and a mode/row-equilibrated diagnostic separately in g0, transition and g1, keeping fixed seed cancellation distinguished. A row-equilibrated numerical experiment would leave the continuum equations and basis unchanged; no such experiment or diagnosed scaling defect is claimed here.

Second, add a coupled flat **vector** RHS-to-JVP inverse control. For one component b^i=v_m, prescribe FM=e_i*Delta(v_m)+(1/3)*partial_i grad(v_m) using the independent Cartesian distance Hessian. Its grad-div term couples Cartesian components and Fourier frequencies, unlike the scalar Laplace inverse above. Recover the original component for bx and a transverse component, initially modes0/1/2/4, with true weighted residual and reconstructed value/gradient errors on several ordinary/core-approach points reported separately. The source/operator controls validate the tensor equations; they do not prove the inverse's conditioning/rounding for this coupled discretization. The present private scalar inverse samples one ordinary point and must not be described as a full spatial error bound.

A common weaker row weight w2=[sin(alpha)sin(beta)]² is a mathematically consistent, separately tagged experiment. It is strictly positive on all Gauss nodes, phi-independent and equal to the analytic polynomial(1-z²)(1-zeta²) in raw grid coordinates. Multiplying both residual and JVP by it leaves all discrete roots, physical boundary conditions and the represented basis unchanged; it introduces no modal q^-r division. Apply it consistently in geometry.weight and the preconditioner row scale. The modal FD block factors are independent of this weight. If J_w=W*J0 and M_w=W*M0, then M_w^-1*J_w=M0^-1*J0, and J_w*M_w^-1=W*(J0*M0^-1)*W^-1. Thus exact-arithmetic preconditioned spectra are unchanged while the residual inner product, line-search/stopping norms and rounding can change substantially. The old A/B distance scaling made w6 behave like r³ near a joint focus and rho³ on an ordinary axis; the new t/eta scaling makes it behave like r⁶/rho⁶. A cubic weight would match that historical distance suppression, whereas w2 deliberately constrains the core more strongly. At N128 the joint-corner w2 is about2.3e-8 versus w6 about1.2e-23. A1e-14 absolute target permits raw residual roughly4.4e-7 under w2 rather than8e8 under w6 at that row. This is a conditioning/stopping experiment, not a PDE or physical-tolerance change; compare unweighted operators through fresh isolated workers and retain the same independent physical gates.

Subsequent read-only inspection of actual saved diagnostics materially qualifies that hypothesis. update_diag divides **every** residual row by its positive weight, so its unscaled_linf is a true global conformal collocation maximum. In converged moderate_c2block N128 its scalar/vector maxima are3.51e-11 and at most9.90e-10; at N152 they are2.52e-10 and at most3.96e-9. Corresponding off-grid bulk physical M RMS values are2.765e-4 and1.744e-4. Large actual raw residuals left unsolved at core collocation nodes are therefore absent in these snapshots. The enormous permitted corner stopping error is a theoretical weakness, not the observed error. A w2 experiment can still diagnose rounding/conditioning sensitivity, but current evidence favors interpolation/source underresolution over a simple early-stopping explanation. No claim that row scaling fixes the binary is made.

Actual shell resolution is more informative than the limiting endpoint spacing. On the smaller moderate core lo/hi=.03/.06,b3, the radial eta1 axis has t(r)=r/(2b+r) and fractional node count(2N/pi)*asin(sqrt(sigma(t))). The angular t0 axis has eta(r)=1-r/b and fractional endpoint count(N/pi)*acos(atanh(eta*tanh(kappa))/kappa). At N80, default .2/2 radial counts inside lo/hi are7.99/11.20 and angular8.91/12.04. Focused .05/3 gives radial15.60/21.42 and angular15.72/19.19. It therefore improves the radial transition spacing by about1.8, but the angular transition count only by about11%, although many more points lie in the already flat inner core. Actual integer Gauss counts round these fractional values with the half-node offset. The angular shell counts at kappa2/2.5/3/4/6 are3.13/3.52/3.47/2.68/1.62; further increasing kappa does not improve this shell. If focused tests identify this bottleneck, independent N_eta refinement is preferable to blindly increasing kappa. These are analytic coordinate-count calculations, not a numerical solve or acceptance result.

Another useful no-solve diagnostic is to prolong a retained converged modal P tensor to a finer grid with the **same map**, then evaluate its native equation samples without solving. The old-node collocation defect compared with the finer-node defect directly exposes off-collocation PDE aliasing independently of Cartesian finite differences. Scalar/source/connection contributions can be separated if needed. Also extend the random-P collocation/sampler comparison to all six Cartesian Hessians using an independent coefficient-based Cartesian chain. Current random-P evidence covers values/first derivatives, while the distance Hessian controls mostly use low-degree/constant P. Passing those does not by itself bound the rounding of a general high-degree retained tensor. Both are proposed private diagnostics and are not counted as completed.

A smaller proposed fallback than domain decomposition recovers original A/B puncture spacing while keeping the modal factors: s=(A+1)/2,zeta=-B,t=s²,eta=2zeta/(1+zeta²),q=s*(1-zeta²)/(1+zeta²), and v=-2(1-s²)q^r P(s,zeta). For r0/1 only, C2 ordinary-axis regularity requires P_s(0)=0 and P_zeta(±1)=0. Frequencies with r>=2 need no endpoint lift: their first extra rho³ harmonic is already C2. Radial and angular maps then have zero endpoint Jacobians and recover O(N^-4) nearest-focus spacing. These zero Jacobians are acceptable only with the stated modal endpoint conditions; blindly changing the grid without the conditions recreates axis cusps.

Only rank1 radial and rank2 angular derivative lifts are necessary for this factored fallback, unlike the rejected rank3/6 V lifts applied to every frequency. Minimal radial null polynomial K=T_N(2s-1) gives p_lift=p-K*p_s(0)/[2*(-1)^(N-1)*N²]. Angular null polynomials K0=T_N(zeta),K1=zeta*T_N(zeta) give p_lift=p+K0*a+K1*b, where [K0'(-1),K1'(-1);K0'(1),K1'(1)]*[a,b]^T=-[p'(-1),p'(1)]^T. Its columns are[(-1)^(N-1)*N²,N²] and[(-1)^N*(N²+1),N²+1], with condition approaching1 after column normalization. They vanish at all original Gauss nodes, retaining primary nodal values; added degrees are N/N+1. Use the direct lifted D2 rather than D1² and exactly the same extension in sampling. Both r1 sine/cosine partners share the extension, so mixed/Fourier derivatives commute. The factors keep v=O(1/r) at infinity. This establishes a consistent C2 interpolation space, **not** elliptic stability: original-grid modal FD chains and first-row regularity still need independent audits, and tiny scalar/coupled-vector Jacobian/inverse controls must precede any production solve. Uniform-angle reflection alone does not prove the additional P_s/P_zeta conditions. No implementation or successful result for this fallback is asserted.

Near-axis sampling for this proposed fallback must evaluate the symbolic Neumann cancellation stably. Forming a rounded endpoint moment p_s and then dividing by s can amplify noise despite an exact-axis clamp. Derivative-basis divided differences, using p_s(s)=p_s(s)-p_s(0), or an explicit stable local polynomial quotient must enter the Cartesian chain before small-coordinate division. The current mapped t/eta sampler's excellent near-axis oracle standard should be preserved; the proposed fallback has not yet met it.

The focused-map library is identified in focused_map_operator_audit.json as SHA e7a88824a3cc8c231ef58af941d86f77350d36895702d3be6ec6251b9c65a65d, exact identifier modal_P_C2prolate_map_v3_r0.050000000000000003_k3. regular_operators_focus05_k3.json now includes scalar and **all three** vector Cartesian components:112rows over N32/64, with largest finest normalized conformal defect1.491e-10. Every finest row passes; some coarse mapped-nonpolynomial rows fail and remain preserved. Seventeen Python controls pass. The separate strict scalar-inverse aggregate remains false as documented above. These successful operator controls do not promote the focused-map binary.

Read-only inspection of the completed moderate-focus05-k3.log shows grids80²/104²/128², nphi28, all terminate with retained Newton line-search failure, weighted residual maxima1.271e-14/2.148e-14/2.894e-14. The near momentum RMS values are3.062e-4/1.135e-3/2.511e-4; bulk1.095e-3/.0101254/.00261068. Charges are unstable, including extrapolated Jz about.5403/6.233/-.1830. These are internally unfinished and physically failed diagnostics, not a resolved convergence sequence. Additional endpoint clustering by itself has not fixed the construction.

The parent's no-solve dense replay provides a stronger localization than the earlier hypothesis. dense_equations_focus05_k3.json restores the same unfinished N80 retained polynomial and evaluates it on104²×28 with no solve. In g1, scaled physical-equivalent My/Mz RMS are4.462e-4/4.422e-4; raw conformal g0 momentum extrema are145–148, compared with roughly2e-9 at the old nodes. dense_equations_phi_focus05_k3.json instead replays80²×56: g1 physical-equivalent momentum component RMS remains about1e-12–5e-12, raw g0 extrema about3e-9. The parameterization/maps/SHA are bound in both reports; no free-source baseline is subtracted. Their interior outputs are scaled **modified equations**, not physical vacuum constraints, while g1 equivalents are physical constraints. This directly demonstrates a large continuous-meridional PDE defect between collocation nodes, independent of the Cartesian FD verifier, without comparable phi aliasing. It does not by itself distinguish source underresolution from poor meridional elliptic conditioning. A retained source's unconverged1.27e-14 weighted maximum cannot explain this many-orders-of-magnitude replay gap.

If row scaling does not alter that defect, a separately declared moderate actual-metric-operator (inner_flatten0) control can help distinguish basis/source resolution from optional-operator conditioning before another major rewrite. The literal Eq27/28 gGamma construction is not metric-compatible and does not inherit the continuum metric-compatible Navier energy proof. Its interior lower-order conditioning is therefore a legitimate test variable. This changes the interior PDE and needs explicit new configuration/evidence; no existing optional-operator outputs should be reinterpreted as that control, and no claim that the published optional attenuation is incorrect is made here.

The newly written private random-polynomial and coupled-vector inverse controls are independently source-audited. cartesian_polynomial builds the Chebyshev recurrence with T0=1 and halves its stored constant only after all recurrences; Cartesian G/D^(r+1), parity-cap, orthonormal Fourier and Nyquist factors match the represented polynomial. It shares the tested AD Jet engine but neither the chart transform nor nodal differentiation matrices. The far0 fixture has no omitted fixed W. random high-degree values/all-six-Hessians checks at N12/24/40 give normalized maxima2.03e-16/8.02e-15/1.75e-13. The exact-delta coupled Navier RHS uses Hessian[d+1][component]/3 plus the Laplacian in the selected bx/by component, with g0,psi1,zero free sources and zero base fields. Sixteen phi nodes resolve the source frequencies throughm6 for inputm4. All24 inverse rows (N8/16/32,m0/1/2/4,bx/by) pass their separately declared relative1e-11,max300Krylov and sampled-field1e-11 bound, using one ordinary and two core-approach points. The largest logged field/gradient error is1.648e-14 at N8,m0,bx; largest N32 is6.122e-15 at m0,by, rather than a universal1e-16 maximum. Most rows are smaller. These checks close the stated manufactured gaps within their tested resolutions; they do not bound an arbitrary N128 binary tensor.

The parent begins a new moderate_focus_wide_actualop case using .05/3 maps, actual metric operators (inner_flatten0), and g radii .2/.8 times the minimum contracted seed Kerr throat. This changes **both** the interior operator and g widths. It can establish a separately supported configuration if physical/charge/axis/covariance/horizon gates pass, but cannot by itself attribute an improvement to either changed ingredient. The .8 seed-throat bound is only an analytic-seed geometric margin until the actual solved binary horizons are found and refined. No accepted result for this new case is asserted here.

The user's enclosure request must be addressed precisely. g and the optional interior operator alter the constraint equations; their support must be inside accepted binary horizons before claiming exterior physical vacuum. The exponential f and F attenuations instead choose free conformal data, have noncompact tails, and cannot literally be wholly enclosed by a finite horizon. Where g1, the actual coupled equations still enforce vacuum even if f/F differ from1. Report their profiles/region residuals, verify physical exterior constraints, and state this interpretation of "modified regions" explicitly rather than claiming all attenuation tails are hidden. A literal requirement to enclose every nonzero f/F tail is not satisfied by the published exponential construction. This distinction is not a waiver of the user's horizon gate or acceptance of the current failed binary.

A cheaper continuous enclosure bound for a represented horizon shape improves the earlier first-derivative Lipschitz estimate. Write R=sum_(l,m)a_lm*Y_lm with orthonormal harmonics and lambda_l=l(l+1). On the unit sphere, integration by parts/Bochner gives integral|Hess_S Y_lm|²=lambda_l*(lambda_l-1). Rotational invariance of the summed addition theorem then gives sum_m|Hess_S Y_lm|²=(2l+1)*lambda_l*(lambda_l-1)/(4pi) pointwise. Cauchy–Schwarz and the triangle inequality yield B2=sum_l||a_l||₂*sqrt[(2l+1)*lambda_l*(lambda_l-1)/(4pi)] as a global Frobenius (hence operator) Hessian bound for R. At its true minimum grad_S R=0; a geodesic delta-cover therefore gives R_min,true >= min_sample R - B2*delta²/2. This quadratic-cover error can be much smaller than the earlier linear Lipschitz bound. For positive-m stored coefficients with real-shape conjugacy, ||a_l||²=|a_l0|²+2sum_(m>0)|a_lm|²; verify the file's actual normalization before use.

For sampled theta rings and uniform phi, an explicit conservative cover is delta=sqrt(theta_cover²+(pi/Nphi)²), with theta_cover=max(theta_first,pi-theta_last,.5*max_ring_gap). It bounds the length of a linear coordinate path because sin(theta)<=1; no pole derivative evaluation is required. The g ball around a seed center is enclosed by the represented star-shaped surface if gmax+|seedcenter-findercenter| < min_sample R-B2*delta²/2, with a numerical margin. This certifies the represented surface only. Strict expansion RMS/max, angular/shape refinement and the physical constraints remain necessary to identify it as the accepted binary horizon. No such enclosure calculation has been run or accepted by the reviewer.

The parent proposes a separate focus-corner regularity diagnostic, not yet implemented. For full r=m, set H=(1-t)^(m+1)*P in v=-4b*G_m*D^(-m-1)*P. Near the plus focus, dx=x-b, r=r_plus and delta=1-eta satisfy t=(r+dx)/(4b)+O(r²),delta=(r-dx)/(2b)+O(r²). The leading nonanalytic G_m*r coefficient is [H_t-2H_eta]/(4b). Thus Cplus[P]=P_t-(m+1)P-2P_eta=0 at(0,1); at the minus focus Cminus[P]=P_t-(m+1)P+2P_eta=0 at(0,-1). These moments are correct necessary analytic-extension conditions for each full-m component in a **flattened g0 core**. Removing them is not sufficient for analyticity. Scalar m0 C2 also requires H_tt-2H_t-4H_etaeta=0 at either focus; its r*dx term is the remaining nonpolynomial degree2 contribution. Mode1 needs the first moment for C2; m>=2's G_m*r is already C2 but is excluded by a genuinely analytic homogeneous core solution.

Minimal proposed radial null lifts Kplus/minus=T_N(z(t))*(1±eta)/2 preserve every Gauss nodal value and ordinary-axis C2. With sN=(-1)^N and p=m+1, their corner-functional matrix is sN*[[-M,1],[1,-M]], M=2N²/lambda+p+1, since T_N(z(0))=sN and partial_t T_N(z(0))=-sN*2N²/lambda. Its condition is(M+1)/(M-1). Physical eta factors are not linear polynomials in the tanh-map raw zeta; the continuous sampler must retain the analytic lift coefficients rather than approximate these factors by the pre-existing angular Chebyshev coefficients. Using raw-zeta linear factors instead replaces the offdiagonal1 by1/eprime and the diagonal p+1 byp+1/eprime, where eprime=eta_zeta(±1). Both tiny matrices are well conditioned for current grids, but this does not establish stability of the resulting PDE interpolation space. Direct lifted D2, mixed-derivative consistency and full scalar/coupled-vector Jacobian/inverse checks remain required before any binary claim.

Scope cautions are essential. Actual Kerr conformal seed metrics are generally only C2 at the puncture (their local expansions include r³ metric terms); an analytic focus extension is not established for inner_flatten0, so high-m analytic conditions must not silently restrict that different curved-core problem. Cappedm>4 require more than these Robin moments: a fully smooth mode needs P~[t*(1-eta²)]^((m-r)/2). A nonzero capped corner P still gives a leading q^r*cos/sin(mphi) nonpolynomial term even if the first Robin condition holds. An initial flattened-core diagnostic can therefore constrain fullm0through4 and retain explicit higher-cap limitations, rather than claiming all focus regularity is fixed.

For the scalar, initially require far_radius0 as well. The harmonic correction in a flattened g0 region is full u=W+v. With F>0, W=(1-F)*(psi_seed-1) starts as r_lab³ times angular dependence at a boosted seed and is generally C2 rather than analytic. Then v must cancel that known nonanalytic W. Homogeneous analytic moments on v would impose the wrong condition: a scalarm2 term can contain G2*r, and a scalarm4 term can contain G4/r, corresponding to an unbounded P_m4~1/r near the excluded focus. An inhomogeneous moment or different scalar correction split would be needed. This is a representation/regularity caveat for far-filter configurations, not a change to the exact W algebra or an explanation of the current far0 failure. Vector b has no W, but the operator/cap qualifications still apply. No corner-lift implementation or acceptance result is counted here.

The old constant-P distance fixture is **not** compatible with these corner moments: Cplus/minus[Pconst]=-(m+1)*Pconst. A newly constrained inverse test must therefore use a different manufactured field rather than expecting the old analytic RHS to recover an inadmissible field. For sigma=t/(lambda+(1-lambda)t), choose P=1+lambda*(m+1)*sigma for fullm>=1; it satisfies both first moments and is exactly low degree in the actual mapped coordinate. For scalarm0 choose P=1+lambda*sigma+lambda*sigma², whose corner Pt1,Ptt2 also satisfy the extra C2 moment. The physical v=-4b*G_m*P/D^(m+1) is C2 at both foci and decays1/r; an artificial manufactured RHS is allowed to be nonzero in the flat core. Compute that RHS with independent Cartesian distance/sigma Jets, not the lifted chart matrices. These fixtures test the constrained interpolation space without a mapped nonpolynomial convergence confound; they are not fully analytic at all higher focus orders and are not vacuum core solutions. A fully analytic rational Cartesian fixture needs its own resolution study. No new fixture has yet been executed by the reviewer.

Checkpoint review frozen on 2026-10-02 for the parent-owned commit. The default native image is SHA 9cbf1108d9060095cade976c27394d8953594f45df97503a90700474f4f822fb. The parent reports the full native/Python suite and refreshed seed gate passing. Read-only inspection of target_seed_controls_current.json confirms this SHA and passing isolated spin95, boost885 and generic combined spin95/boost885 records, with historical proof preserved separately; these are analytic-seed controls, not solved-binary acceptance. Read-only inspection of collocation_map_api_equivalence.json confirms separate controlled processes, verified loaded paths and exact equality of the retained off-axis fields, equation residual and JVP between c918 and 9cb. This is the read-only map-query/basis-identifier migration, distinct from the earlier derivative-rounding comparison and the invalidated same-process dyld comparison. Its axis finite-difference constraints remain finite but physically unacceptable; equality of two builds does not improve those constraints.

The additional no-solve meridional replays retain the exact same focused N80 polynomial and the same bound SHA/maps. On 80×160×28, g1 physical-equivalent My/Mz RMS are 1.2467135e-4/1.2934113e-4; on 160×80×28 they are 2.9810327e-5/3.0704006e-5. Phi-only 80²×56 remains below about 5e-12. Thus both meridional directions expose a substantial continuous PDE defect, with the angular replay larger in this diagnostic. These are evaluations of the same interpolant, not separately solved refinement sequences or acceptance. The wide-g actual-operator N56 test terminates with a Krylov failure; its N80 test is deliberately interrupted. The matched wide-g flattened N80/N104 outputs also fail physical acceptance. No corner lift, row-weight change or accepted generic binary is introduced at this checkpoint. An asymmetric polar refinement is only planned as a separately declared experiment.

The private constant-anisotropic-metric test is independently source-audited. With unit n=(.5,.8,.3)/sqrt(.98), h=I+epsilon*n*n^T and h_inverse=I-epsilon/(1+epsilon)*n*n^T are exact Sherman–Morrison pairs. Constant background Jets, psi1, zero connection/free sources, zero base fields and complete cache replacement isolate the differential operator; only the original row weights are retained. The scalar RHS is h^{jk}v_jk. For a single nonzero contravariant vector component b^c=v, the vector RHS is delta^i_c*h^{jk}v_jk+(1/3)h^{ij}v_jc. This matches both the independent distance-Jet fixture and the cached L/divergence implementation, including the Hessian index order. The fixture shares the separately tested Jet engine but does not reuse the chart Hessian or spectral differentiation matrices.

anisotropic_inverse_default.json binds the original test source SHA 50d3d5c242c969aceb8e55c94636288e5857f5c095fff34bce16264d35fd3ae3 and executable SHA c0ad1bb8a33b6f9e3014b8682682a461accd8d9177b19444f008862d08bb7e41 to the 9cb image. The declared audit has 72 rows: epsilon .01/.4, N8/16/32, modes0/1/2/4 and scalar/bx/by, restart32, max300, relative tolerance1e-11 and sampled field/gradient bound1e-11. Its aggregate remains **false** with seven failures, all epsilon .4, N32 vector rows. In particular m0 bx/by finish at relative residual1.27271e-4/1.60582e-4 and field/gradient error8.89912e-5/1.23788e-5. All epsilon .01 rows and the epsilon .4 scalar rows pass. The original source/executable/log remain archived.

anisotropic_inverse_restart128.json is a separate, passing diagnostic of only those first two failed rows, with identical geometry/RHS/PDE/preconditioner/tolerances, restart128 and max1024. Its changed test source SHA is f6376af06e6d510b021bdee632f718fb6a83273149497a96fc2a9113e9ed0c13. bx/by converge in82/85 iterations with independently recomputed true relative residual9.23627e-12/9.89474e-12 and field/gradient error1.46506e-12/5.56348e-15. This demonstrates restart stagnation for these two RHSs; it does not prove all original failed rows pass, full matrix rank, or a binary accuracy improvement. No original failed flag is replaced.

For constant h on a compatible whole-space decaying/regular domain, define Pi_long=grad_h*Delta_h^{-1}*div. Commuting derivatives imply Pi_long²=Pi_long and L_h=Delta_h*(I+Pi_long/3), hence L_h^{-1}=(I-Pi_long/4)*Delta_h^{-1}. This gives a principled coupled Helmholtz preconditioner proposal using scalar anisotropic Poisson inverses plus divergence/gradient operations. In the h inner product the continuum longitudinal/transverse symbol ratio is4/3; epsilon .4 gives metric eigenvalue ratio1.4, so comparison with an exact isotropic scalar inverse has a modest continuum bound of(4/3)*1.4. That bound does not apply to the weighted modal coordinates, approximate FD inverse or finite capped interpolation domain. Discrete derivative/inverse commutation and the regularity/boundary domain must be checked independently before calling a finite implementation exact. The passing larger-restart diagnostic weakens a claimed discrete-nullspace explanation for its tested RHSs; it neither establishes a stable global discretization nor fixes the observed meridional binary defect. No Helmholtz implementation or numerical result is counted here.

**Proposed stable scalar-source grouping, review revision of 2026-10-02.** The derivation below is complete as an algebraic proposal. It has not been implemented, numerically validated or accepted. No source change, solver job or new physical acceptance is attributed to this review. The first proposed implementation retains the original Ψ/F split, the fixed W contribution inside u, and the existing correction operators. The advanced curvature and W reorganizations described later are optional future measures.

The motivation is the retained failed infinity-equilibration control, not a claim that the original physical seed fields fail. In far_source_floor_infinity_equilibrated.json, image SHA 40672d36b1a1d1d3deda9a9fb38ea4f243da05c4dfde4dea30856e651aa58dd4 uses sin6/(1−t)^6 and declares a weighted exact-seed far-source bound of1e-14. At 80×160×28, boost885 and combined_generic reach scalar floors1.8822e-13/1.8412e-13, while their raw physical-equivalent scalar extrema are about3.5e-18. That new-norm gate is false and no binary solve in it is accepted. The original norm and binary gates remain in force. The proposed regrouping must undergo new-SHA controls before any renewed equilibration experiment.

All sums below run over active seeds; there is no branch selecting an isolated seed. Write p_s=ψ_s, h_s for its conformal metric, A_s for its covariant conformal trace-free tensor and K_s for its mean curvature. Set w_s=f_s F_s and use the actual background fields

$$
h=\delta+\sum_s w_s(h_s-\delta),\qquad
K=\sum_s w_s K_s,\qquad
S=\sum_s A_s,\qquad
\tau=\operatorname{tr}_h S,\qquad
M=S-\frac13h\tau .
$$

Indices of a contraction with subscript h are raised twice with h inverse; Q_s=|A_s|²_{h_s}. The seed identity, with the existing curvature and extrinsic-curvature signs, is

$$
\Delta_s p_s=
\frac18 p_s R_s+\frac1{12}p_s^5K_s^2
-\frac1{8p_s^7}Q_s .
$$

It follows from the exact Kerr seed vacuum Hamiltonian constraint in either supported conformal choice. The construction must retain the actual metric's Levi-Civita connection for every seed and source identity. The optional opGamma=g Gamma is not used in this identity or in the source curvature.

For the first implementation define Z=Ψ=1+sum_s F_s(p_s−1), ψ=Z+u and L=longitudinal_h(b). The current u includes W in all ten scalar field slots; that convention remains unchanged. The complete background intercept to cache is

$$
C_Z=\Delta_h Z-\frac18 ZR_h-\frac1{12}Z^5K^2
+\frac{|M|_h^2}{8Z^7}.
$$

Replacing only lapPsi by the seed identity does not solve the problem: the runtime kernel would still subtract its large curvature and seed terms. The complete C_Z is the object to assemble with the following differences.

For each seed form d_s=h−h_s directly from the filter contributions, as in the existing momentum source:

$$
d_s=(w_s-1)(h_s-\delta)+\sum_{t\ne s}w_t(h_t-\delta),
\qquad
e_s=h^{-1}-h_s^{-1}=-h^{-1}d_s h_s^{-1}.
$$

Form D_s=Gamma_h−Gamma_s stably from d_s, rather than subtracting two rounded connections:

$$
(D_s)^i{}_{jk}=\frac12h^{il}
\left[(\nabla_s)_j(d_s)_{kl}
+(\nabla_s)_k(d_s)_{jl}
-(\nabla_s)_l(d_s)_{jk}\right].
$$

The exact Laplacian difference on p_s is

$$
E_s=(\Delta_h-\Delta_s)p_s
=e_s^{ij}\left[p_{s,ij}-(\Gamma_s)^k{}_{ij}p_{s,k}\right]
-h^{ij}(D_s)^k{}_{ij}p_{s,k}.
$$

All F derivative terms in Δ_h Ψ are retained:

$$
T_F=\sum_s\left[
2h^{ij}F_{s,i}p_{s,j}+(p_s-1)\Delta_h F_s
\right],
\qquad
\Delta_h Z=\sum_sF_s\left[\Delta_s p_s+E_s\right]+T_F .
$$

There is no f multiplier in Z or T_F. Derivatives of f enter h, D_s and the actual curvature through w_s; derivatives of F enter those fields and T_F. Multiplying the source brackets by g does not introduce derivatives of g.

The grouped intercept is

$$
C_Z=T_F+\sum_sF_s E_s+\frac18 C_R+\frac1{12}C_K+\frac18 C_A,
$$

with exact blocks

$$
\begin{aligned}
C_R&=\sum_sF_s p_s(R_s-R_h)+\left(\sum_sF_s-1\right)R_h,\\
C_K&=\sum_sK_s^2\left[F_sp_s^5-w_s^2Z^5\right]
-2Z^5\sum_{s<t}w_sw_tK_sK_t .
\end{aligned}
$$

For the tensor norm, construct its metric change without subtracting two norms:

$$
dQ_s=\left(e_s^{ik}h^{jl}+h_s^{ik}e_s^{jl}\right)
(A_s)_{ij}(A_s)_{kl}.
$$

Since each seed is analytically trace-free in its own metric, τ=sum_s e_s:A_s is the stable trace already used in the momentum source. Then

$$
C_A=
\frac{\sum_s dQ_s+2\sum_{s<t}\langle A_s,A_t\rangle_h-\tau^2/3}{Z^7}
+\sum_sQ_s\left[Z^{-7}-F_sp_s^{-7}\right].
$$

The −τ²/3 projection term is exact analytically. For closer consistency with a supplied rounded projection tensor, set C=−hτ/3 and replace the cross/projection numerator by

$$
2\sum_{s<t}\langle A_s,A_t\rangle_h
+2\langle S,C\rangle_h+|C|_h^2.
$$

This follows from M=S+C and avoids assuming that a freshly computed floating-point trace of S is bitwise τ. Use common contraction helpers/order in the cached and point kernels. Raw trace, inverse and projection discrepancies should be recorded as arithmetic errors; they must not be hidden by weakening a physical gate.

The curvature difference is also formed covariantly. With D_s stored as Jets,

$$
(\nabla_s)_p D^i{}_{jk}
=D^i{}_{jk,p}
+(\Gamma_s)^i{}_{pl}D^l{}_{jk}
-(\Gamma_s)^l{}_{pj}D^i{}_{lk}
-(\Gamma_s)^l{}_{pk}D^i{}_{jl}.
$$

The Ricci convention matches the current curvature helper:

$$
\begin{aligned}
\delta\operatorname{Ric}_{ij}
&=(\nabla_s)_kD^k{}_{ij}
-(\nabla_s)_jD^k{}_{ik}
+D^k{}_{ij}D^l{}_{kl}
-D^l{}_{ik}D^k{}_{jl},\\
R_h-R_s&=e_s^{ij}(\operatorname{Ric}_s)_{ij}
+h^{ij}\delta\operatorname{Ric}_{ij}.
\end{aligned}
$$

Thus C_R uses the negative of this δR for R_s−R_h. Only first derivatives of D_s are required. Existing metric AD2 and seed connection first derivatives suffice: diff(metric) followed by a connection-difference Jet supplies the needed slots. Its absent or fabricated higher derivatives must not be used. The existing seed_sum_source already forms seed inverse, connection, d_s and e_s, but its current DC stores values only; scalar curvature needs a Jet D_s with first derivatives. Reusing those temporaries is more practical than repeating all seed geometry work.

Power differences must be factored. In particular,

$$
\begin{aligned}
F_sp_s^5-w_s^2Z^5
&=(F_s-w_s^2)p_s^5
-w_s^2(Z-p_s)\sum_{k=0}^4Z^{4-k}p_s^k,\\
F_s-w_s^2
&=F_s\left[(1-F_s)+F_s(1-f_s)(1+f_s)\right],\\
Z-p_s
&=(F_s-1)(p_s-1)+\sum_{t\ne s}F_t(p_t-1).
\end{aligned}
$$

Also Z^−7−F_s p_s^−7=(1−F_s)p_s^−7+(Z^−7−p_s^−7). The last difference can use the exact factored rational polynomial

$$
Z^{-7}-p_s^{-7}
=-\frac{(Z-p_s)\sum_{k=0}^6 Z^{6-k}p_s^k}{Z^7p_s^7},
$$

or expm1(−7 log1p((Z−p_s)/p_s))/p_s^7 with positive p_s,Z and a suitable overflow-safe ratio evaluation. This is a proposal for stable arithmetic, not a tested choice of implementation. Near-unity filter differences may likewise benefit from expm1, but changing their finite-precision Jet evaluation requires explicit field/derivative comparison; their mathematical definitions do not change.

With C_Z cached, the runtime nonlinear correction increment is

$$
\begin{aligned}
I(u,L)
={}&-\frac18uR_h
-\frac1{12}K^2\left[(Z+u)^5-Z^5\right]\\
&+\frac{2\langle M,L\rangle_h+|L|_h^2}{8\psi^7}
+\frac{|M|_h^2}{8}\left[\psi^{-7}-Z^{-7}\right],
\qquad \psi=Z+u,\\
F_H={}&\Delta_{\mathrm{op}}u+g\left[C_Z+I(u,L)\right].
\end{aligned}
$$

Use (Z+u)^5−Z^5=u sum_(k=0)^4 (Z+u)^(4−k)Z^k and the corresponding inverse-power difference rather than two power evaluations followed by subtraction. The analytic JVP remains the existing one:

$$
\delta F_H=\Delta_{\mathrm{op}}\delta u
+g\left[
\left(-\frac18R_h-\frac5{12}\psi^4K^2
-\frac7{8\psi^8}|A|_h^2\right)\delta u
+\frac{\langle A,L(\delta b)\rangle_h}{4\psi^7}
\right],\qquad A=M+L.
$$

The intercept is independent of the correction. Physical ψ, h, K, M, sampling and the optional correction operators retain their existing definitions. In the isolated unfiltered limit, d_s,e_s,D_s,T_F and all coefficient differences vanish, so the complete intercept cancels as a limit of the generic formula. No isolated-seed if branch is needed. This is equality of the continuum PDE; floating-point residuals may change and bitwise old/new residual equivalence must not be claimed.

The required raw-oracle proof explicitly tracks rounded seed identities. Define

$$
\epsilon_s=\Delta_s p_s-\frac18p_sR_s-\frac1{12}p_s^5K_s^2
+\frac{Q_s}{8p_s^7}.
$$

In exact algebra, the original direct residual minus the grouped residual is g sum_s F_s ε_s. Floating-point metric/inverse, trace/projection and reassociation errors add to that measured difference. Retain ε_s and verify the discrepancy against this prediction with an independently stated arithmetic bound. A smaller grouped residual alone is not proof: the seed identity sets a known analytic zero, whereas its rounded raw geometric evaluation has a nonzero floor.

Required controls, all still proposed, are:

- Preserve an untouched raw path that computes actual Δ_h Ψ, Ricci/scalar curvature, total projected A and the nonlinear source directly, without the seed vacuum substitution. Verify connection and Ricci differences against that path before testing the complete grouping.
- Use generic unequal seeds with arbitrary spin/boost directions, both conformal choices, nonzero scalar/vector correction Jets, f/F transition derivatives, g0/transition/g1 and both correction-operator choices. Include ordinary, far, inner-sheet and regular throat points within their valid domains. Compare cached and point kernels with common conventions.
- Supply a private nonvacuum seed perturbation, such as a controlled change to scalar p_s Jets while retaining a trace-free A_s, and restore the measured ε_s terms in the grouped formula. This negative control must recover the raw source. A mismatch would expose an invalid general algebra or an oracle that merely imposes vacuum.
- At g1 independently evaluate raw physical geometry and its Hamiltonian constraint, using H_phys=−8ψ^−5 F_H in exact arithmetic. Retain the Cartesian finite-difference seed controls, field/gradient comparisons, correction JVP finite differences and raw geometric error floors. At g<1 distinguish modified equations from physical vacuum constraints.
- Bind all evidence to the new library/source SHA. Only after the raw-oracle proof passes, rerun the declared infinity far-source gate and the analytic flat/coupled Navier inverse controls in that new norm. Preserve the failed earlier equilibration artifact. No cross-SHA binary acceptance, charge acceptance or relaxed tolerance follows from this derivation.

Two optional future reorganizations may reduce remaining cancellation, but neither is required for the first Ψ/F implementation. First, per-seed R differences can still subtract leading binary far terms. Rewrite

$$
C_R=\sum_sF_s(p_s-1)(R_s-R_h)+\left(\sum_sF_sR_s-R_h\right).
$$

For a Cartesian metric h=δ+q define its flat-linear scalar curvature L_flat(q)=partial_i partial_j q_ij−Delta_flat tr(q), and R(h)=L_flat(q)+N(h). Build Gamma_lin=(1/2)δ inverse times the derivative permutations of q and Gamma_nl=(1/2)(h inverse−δ) times those permutations, with h inverse−δ=−h inverse q formed stably. N comprises (h inverse−δ):Ric_lin plus h inverse contracted with the derivatives of Gamma_nl and the full Gamma*Gamma terms. This separates additive O(r^−3) far terms from nonlinear O(r^−4) terms. With q_s=h_s−δ,

$$
\begin{aligned}
\sum_sF_sR_s-R_h
&=\sum_s\left[(F_s-w_s)L_{\rm flat}(q_s)
-T_{\rm flat}(w_s,q_s)\right]
+\sum_sF_sN(h_s)-N(h),\\
T_{\rm flat}(w,q)
&=2w_{,i}\partial_jq_{ij}+w_{,ij}q_{ij}
-2w_{,i}\partial_i\operatorname{tr}(q)
-(\Delta_{\rm flat}w)\operatorname{tr}(q).
\end{aligned}
$$

For f=F=1, the additive flat-linear terms cancel before floating-point evaluation. This advanced decomposition needs its own direct-curvature proof and seed/binary far-floor measurements; no benefit is asserted as already demonstrated.

Second, the exact fixed-W reparameterization permits Z=Ψ0=1+sum_s(p_s−1), u=W+v and W=sum_s(1−F_s)(p_s−1). The generic intercept formulas then use scalar weights a_s=1 instead of F_s, while metric/mean-curvature weights remain w_s=f_sF_s and the F-product term T is absent. In that split

$$
F_H=\Delta_{\rm op}v+
\left(\Delta_{\rm op}-g\Delta_h\right)W
+g\left[C_{\Psi_0}+I(v,L)\right].
$$

The W defect must be evaluated from coefficient differences:

$$
\left(\operatorname{opinv}^{ij}-g h^{ij}\right)W_{,ij}
-\left[
\operatorname{opinv}^{ij}(\operatorname{opGamma})^k{}_{ij}
-g h^{ij}(\Gamma_h)^k{}_{ij}
\right]W_{,k}.
$$

It vanishes at g1; for inner_flatten0 it is (1−g)Δ_h W, whereas the literal inner_flatten1 connection is opGamma=g Gamma_h and needs its own contraction. This can avoid a later cancellation between Δ_op W and g Δ_h Ψ_F, but requires consistent correction-field conventions throughout cached, point and sample paths. It is deliberately deferred from the first implementation. The original fixed-W convention remains the basis of the required initial raw-oracle proof.
