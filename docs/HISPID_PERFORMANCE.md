# HiSpID serial performance and integration handoff

The solver and the AthenaK importer are implemented in isolated branches.
Exact isolated chi=.95 and v=.885 controls, separately and combined, pass
independent physical checks and AthenaK initial-time horizon controls.
The current regular-basis moderate binary passes preliminary physical and
refined charge, solved coordinate-covariance and direct AthenaK horizon/enclosure
gates. Stronger accuracy and revised high-parameter binaries remain pending. Its high-spin/boost range is
unvalidated; do not use these binary outputs as validated production data.

## Measured improvements

- At identical moderate80²×28 free data, exact block elimination of the
  five-point modal FD preconditioner reduced the solve from296.95s/896
  Krylov iterations to31.33s/64 iterations:9.48x. The full residual and JVP
  remain pseudospectral; this is a preconditioner improvement. Physical
  failure flags are retained (`modal_block_equivalence.json`, `results.json`).
- ADM energy now uses the independently tested analytic physical metric
  gradient. At matched16×32 sphere nodes and R40/200/1000, integration time
  fell13.657s→1.120s (12.20x). Energy differences are<2.7e-11; P/J are
  bit-identical when grid points match (`analytic_charge_method_comparison.json`).
- Aligning the sphere polar axis with the prolate axis removes the severe
  azimuthal alias of the retained meridional polynomial. On the same failed
  80×160×28 binary at R200,256×64→256×128 changes EPJ by6.19e-14;
  256×128→384×128 changes it by6.35e-10. Independent Cartesian-FD fluxes
  and centered tensor covariance pass1e-8 at R40/200/1000; maximum
  discrepancies are5.52e-11 and1.61e-10. This validates integration, not
  the binary geometry (`aligned_charge_fd_comparison.json`,
  `charge_quadrature_polar2_aligned.json`).
- The allocation guard now counts the current five-point stencil instead
  of the obsolete76-point stencil. The default2048MiB and explicit8192MiB
  maximum are unchanged. The new estimate also retains dense radial block
  factors, their shared Fourier/vector groups, and the GMRES basis.
- Centered aligned sphere rings now reuse P/Pt/Peta across phi. Matched
  64×64 quadrature at R40/200/1000 on128×256×28 data falls22.081s→1.046s:
  another21.12x, with maximum EPJ difference5.52e-13. Ordinary sample,
  residual and JVP paths remain bit-identical. Phase-shifted all-mode
  first-gradient, translated/rotated-frame and changed-unknown controls
  pass; independent centered/off-center FD fluxes agree within5.82e-11.
  Off-center spheres retain ordinary interpolation.

A subsequent exact serial reuse removes repeated azimuthal averaging per
Fourier row and moves the vector LU/transfer bank between Newton steps.
Scalar potential/factors and every row scale still rebuild. On a matched
fresh rotated80×160×28 solve, time falls59.97s→55.44s (1.082x in one paired
measurement), with identical4Newton/45Krylov counts and bit-identical final
unknowns, physical fields and residuals. Peak RSS is recorded per worker.
The rectangular private control proves bitwise matrix/factor/random-RHS
equivalence, changing scalar potential, first/second factor counts7+7/7+0,
malformed-cache rejection and direct modal inverse error5.9e-18. A scaled
Fourier round-trip fixture failed identically in cached/uncached paths at
6.5e-10; that first diagnostic remains retained. Separate-process field
checks on allthree saved grids and80×160 residual/JVP checks are bitwise.
All native and17Python controls pass. These source changes do not relabel
the original binary acceptance records.

On the current80×160/104×208/128×256×28 moderate sequence, all near/bulk
physical H/M RMS values improve monotonically. The finest nearH/M are
1.95e-6/6.02e-6 and bulkH/M1.25e-7/8.01e-5; the209.0s solve uses39Krylov
iterations on one thread. Each grid's ADM extraction uses polar orders
2Npolar/3Npolar and phi64/128 at R100/200/400. All satisfy the declared
1e-7 integration bound. The two finest extrapolated EPJ vectors differ by
2.33e-7. `polar_sequence_refined_charges.json` records this preliminary
qualification and its explicit failed stronger/full-binary flags.

Fresh rotated/translated solves of the same sequence pass the declared
coordinate gate. Finest metric/K rotation errors are3.92e-11/2.84e-6;
translation errors are<3.5e-13. Centered charge error is3.40e-10 and fixed
global-origin charge error1.99e-4 versus the preliminary.005 bound. The
fixed-origin extraction still uses coarse12×24 integration and needs
separate angular refinement for stronger claims. The whole covariance
process peaks at10.32GB over sequential contexts (one CPU thread), while
the original finest solve alone peaks at7.83GB. Context budgets do not
bound allocator retention or whole-process RSS.

Jobs and builds use one CPU thread. No OpenMP/MPI/GPU implementation or
scaling claim is made. CPU time and peak memory must justify that next
step. Perlmutter SSH access is available; no allocation has been made.

## Review and replay

Native branch `codex/hispid` lives in the sibling `TwoPuncturesC` worktree,
based on68287742f4920f4ea39b7dac1571c81eefe2ff8f. The BY sources are untouched.
`docs/HISPID.md` contains build, ABI, conventions and replay instructions;
`docs/VALIDATION.md` and retained JSON give configurations and failed gates.
Always supply an absolute native-library path. A checkpoint is bound to its
source SHA, basis identifier, maps and full grid. Normal replay rejects a
mismatch. Explicit separate-process bitwise witnesses may authorize an
API-only initial guess or read-only integration check; fresh solves and
acceptance checks are still required. They never relabel a failed binary.

AthenaK branch `codex/hispid-pgen` lives in the sibling `AthenaK` worktree,
based exactly on PR790 head22baa243970fa1880b2bbc48e88a590069d55e47 from
`HengruiZhu99/athenak:project/z4c_overhaul`. Its `z4c/hispid` pgen loads
portable physical gamma/K data, fills active and ghost cells, and checks
ADM/Z4c round trips. `docs/hispid.md` gives its build/replay commands;
`docs/hispid-current-controls.json` binds the current executable and exact
seed tests to native9cbf1108… source provenance. Direct native-geometry
finder controls and independently refined mesh-constraint checks are
reported separately. These controls take zero evolution steps.

The checked128×256×28 moderate checkpoint has now been exported and tested
with the separate exact-source AthenaK consumer. Both component surfaces pass
expansion RMS1e-7 at lmax16,ntheta32/48. Fixed-order area changes are<7e-12
relative. Real-harmonic coefficient bounds certify complete inner g/operator
balls on the retained surfaces, with margins.176225M/.115817M after an
empirical refinement allowance. Upper-radius bounds certify distinct
components. The mesh import round-trip error is4.1e-16; no evolution steps
are taken. See validation/polar_sequence_horizons.json and the sibling
AthenaK docs/hispid-moderate-binary.json for source/consumer fingerprints.

This direct-geometry finder check is separate from mesh-resolved finder
accuracy. The noncompact f/F attenuation tails require exterior constraint
checks and cannot be described as entirely inside a finite horizon. Reported
coordinate spin is not a generic approximate-Killing-vector spin. Stronger
physical accuracy and revised high-spin/boost binaries remain unvalidated.

No merges, pushes, shared installations, main-project branch changes or
production evolutions have been performed. GPU work remains last priority.

The revised rest-spin chi=.95 binary now has three successful serial Newton
solves. Its128×256×24 solve takes219.32s with4Newton/54Krylov iterations;
the whole solve/verifier process peaks at6.74GB. Physical constraints still
fail, so timing is a diagnostic workload rather than an accepted high-spin
benchmark. The exact126300dc producer exports the labeled checkpoint and
AthenaK imports it with ADM/Z4c round-trip error4.11e-16.

Both coarse horizons are measured, with each Christodoulou mass.5006533,
area8.3135384 and coordinate chi.9475223. A measured-radius initial guess
reduces the matched coarse search213.69s→44.39s (4.81x), with relative area
change9.4e-13; this changes the initial guess only. The complete bounded
measurement/failed-fine-control process takes454.97s and peaks234.11MB on
one CPU. Fine lmax16 expansion stalls above the strict threshold, with
angular variance near2.32e-7. Fixed-order quadrature/higher-order checks
remain pending. See validation/spin95_local_horizons.json. No new parallel
or accelerator work is justified by these timings alone.
