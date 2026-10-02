# Equation system and linear backend selection

The physical equation system and the Krylov algorithm are independent choices.
Both HiSpID (four coupled curved equations) and Bowen–York (one Hamiltonian
constraint with analytic momentum data) support GMRES and BiCGStab. Existing
entry points keep HiSpID/GMRES and BY/BiCGStab as defaults.

`src/PunctureKrylov.c` contains the only production GMRES and BiCGStab iteration
implementations. Its matrix-free callbacks apply a fixed spectral Jacobian and
right preconditioner during each linear solve. The system adapters own the
physical residual, Newton forcing and damping, preconditioner construction and
factor reuse. BY solves `J dv = F` and subtracts `dv`; Hi solves `J step = -F`
and adds the damped step. No PDE, map, attenuation or source term changes.

GMRES retains lazily allocated Arnoldi/preconditioned columns, two-pass modified
Gram–Schmidt, restart and explicit true-residual checks. BiCGStab keeps eight
vectors and checks its true residual on new paths. There is no automatic switch
to another algorithm after failure. Status distinguishes iteration/true-residual
limit, breakdown, invalid inputs, callback failure, nonfinite output and shared
workspace allocation failure. This does not promise exhaustive allocation
handling in the pre-existing system adapters.

BY's original out-of-line `norm2` and `scalarproduct` reductions are supplied as
optional BiCGStab callbacks. GMRES uses its original direct reductions,
keeping Hi's rounding unchanged. This preserves rounding barriers: the Apple compiler fused
square/product accumulation when the same loops were extracted and inlined,
which initially broke bitwise preservation. Hi's original reduction arithmetic
is retained. The default BY recurrence-only stopping policy, absolute breakdown
thresholds and update behavior are preserved; new GMRES, modal or fixed relative
paths fail closed after a failed linear solve.

## Native interfaces

Build serially with `make -j1 all hispid`. Existing Hi config and diagnostic
struct layouts are unchanged. Use the size-tagged options API:

```c
#include "HiSpID.h"
#include "PunctureKrylov.h"
HiSpID_SolveOptions options;
HiSpID_default_solve_options(&options);
options.krylov = PK_BICGSTAB;       /* or PK_GMRES */
options.linear_rtol = 1e-3;         /* 0 keeps native adaptive forcing */
int status = HiSpID_solve_with_options(context, &options);
```

For BY, set options before `TwoPunctures_make_initial_data()`:

```c
TwoPunctures_params_set_Int("TP_krylov_solver", PK_GMRES);
TwoPunctures_params_set_Int("TP_krylov_maxit", 2000);
TwoPunctures_params_set_Int("TP_krylov_restart", 64);
TwoPunctures_params_set_Int("TP_preconditioner", 1);  /* 0 lines, 1 modal */
TwoPunctures_params_set_Int("TP_linear_relative", 1);
TwoPunctures_params_set_Real("TP_linear_rtol", 1e-3);
```

Defaults are BiCGStab, cap100, restart64, lines, and the original absolute
inner target. GMRES can also use the inherited line preconditioner; backend
selection does not select a preconditioner. `TwoPunctures_params_reset()` frees
an unowned table after failed setup and refuses to touch a live solve.

## Python selection

Add `python` to `PYTHONPATH`; all library paths must be absolute. The common
factory dispatches the equation system while `solve` selects the algorithm:

```python
from punctures import Backend
backend = Backend("hispid", "/absolute/build-hispid/libHiSpID.so")
config = backend.config()
with backend.create(config) as data:
    result = data.solve(krylov="bicgstab", linear_rtol=1e-3)
```

```python
backend = Backend("bowen_york", "/absolute/lib/libTwoPunctures.so")
config = backend.config()
config["real"].update(par_P_plus2=.02, par_S_plus3=.03)
with backend.create(config) as data:
    result = data.solve(krylov="gmres", linear_rtol=1e-3,
                        preconditioner="modal", max_krylov=2000, restart=64)
```

Each system keeps its own physical input semantics: Hi seed velocity/rest spin
are distinct from BY lab momentum/angular momentum. BY config has typed
`real`/`integer` dictionaries, and its factory chooses spectral sampling.
Unknown parameter names and invalid Python tolerances are rejected before
native mutation. `None` retains native forcing; an explicit zero is invalid in
Python. `resolved_options` records choices. Archived libraries reject requested
features they lack. Both adapters verify the loaded library image and digest.

BY remains process-global and serial. The Python adapter permits one live
context per process, prohibits a second solve in that context, releases both
data and parameters on close, and checks outer convergence and target masses
rather than interpreting native diagnostic status alone as convergence. For
parallel independent solves use fresh processes. Do not mix raw BY parameter
mutations with a live Python-owned context.

## Evidence and limits

The reproducible driver is `validation/benchmark_krylov_matrix.py`; criteria
were recorded in `validation/krylov_matrix_acceptance.json`. Results retain
fresh-process timings, counts, true linear histories, independent Cartesian
constraint checks, hashes and complete states. The moderate binary uses
40×80×16, one CPU thread, zero guess, cubic residual rows, outer maximum1e-12,
true RHS-relative L2 target1e-3, cap2000 and restart64, with identical modal
preconditioners within each system. Two interleaved repetitions are required.
This is a computational comparison, not equal physical accuracy between BY
and HiSpID.

All four paths meet those stopping rules. BY cross-backend data passes the
core equivalence gate. Hi BiCGStab agrees closely in sampled physical fields,
charges and finite-difference constraints, but fails the declared full raw-P
and coefficient equivalence bound. Its largest raw-P discrepancy lies in m6
near the joint axis/focus, where the physical reconstruction multiplier is
about3.8e-17. This is evidence of weak representation observability, not proof
of an exact null space. A tighter linear tolerance alone did not resolve the
small native-weighted fixture's discrepancy. Hi BiCGStab stays opt-in; do not
claim interchangeable complete saved Hi data or new physical binary acceptance.
The failed v4 extraction experiment and small fixture logs are retained.

`make -j1 test-hispid` includes independent nonsymmetric matrix controls,
callback/nonfinite/breakdown failures, selector/lifetime guards and existing
continuous operator/preconditioner controls. BY original small fixtures,
including the target-mass loop and BL limit, run through
`validation/verify_krylov_backends.py` using both backends. Benchmark state
capture and physical checks occur after solve timing/peak-RSS capture.
