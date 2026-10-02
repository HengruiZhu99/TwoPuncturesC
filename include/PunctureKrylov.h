#ifndef PUNCTURE_KRYLOV_H
#define PUNCTURE_KRYLOV_H
#ifdef __cplusplus
extern "C" {
#endif
/* Equation-independent, serial, matrix-free right-preconditioned solvers.
 * apply/precondition return0 on success. Callbacks may reuse private scratch,
 * but must not retain the supplied vectors; input/output do not alias.
 * A and M remain fixed for the duration of one solve. */
enum { PK_GMRES=0, PK_BICGSTAB=1 };
enum { PK_SUCCESS=0, PK_LIMIT=1, PK_BREAKDOWN=2, PK_CALLBACK=3,
       PK_INVALID=4, PK_ALLOCATION=5, PK_NONFINITE=6 };
typedef int (*PK_Apply)(void *,const double *,double *);
typedef void (*PK_Monitor)(void *,int,double,double,double,double);
typedef struct {
  int method,max_iterations,restart,verify_true_residual;
  int legacy_bicgstab,eager_gmres_basis; /* regression/test policies */
  double absolute_tolerance;
  /* Optional inherited BiCGStab reductions preserve rounding barriers.
   * GMRES retains its original reductions (these fields are ignored). */
  double (*dot)(const double *,const double *,int);
  double (*norm)(const double *,int);
} PK_Options;
typedef struct {
  int status,iterations,operator_calls,preconditioner_calls;
  double recurrence_residual,true_residual,relative_residual;
} PK_Result;
/* x is an initial guess for BiCGStab; GMRES retains its zero-start contract.
 * The legacy policy preserves inherited absolute1e-50 breakdown checks and
 * recurrence arithmetic. All new paths verify the true residual explicitly.
 * Result counters include initial/true-residual operator calls. */
int PK_solve(int n,const double *rhs,double *x,const PK_Options *,
             PK_Apply operator_action,PK_Apply precondition,void *context,
             PK_Monitor monitor,void *monitor_context,PK_Result *);
const char *PK_status_string(int status);
#ifdef __cplusplus
}
#endif
#endif
