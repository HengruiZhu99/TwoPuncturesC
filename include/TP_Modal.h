#ifndef TP_MODAL_H
#define TP_MODAL_H
#include <stddef.h>
/* Fourier/azimuthally averaged block-tridiagonal inverse of the existing
 * finite-difference BY Jacobian. It is a preconditioner only: the spectral
 * residual and Jacobian-vector product are never replaced. Scalar nvar=1,
 * even nphi>=4. Factors/workspaces belong to one fixed Newton Jacobian.
 * Serial/non-reentrant, like the process-global BY solver.
 */
#ifdef __cplusplus
extern "C" {
#endif
typedef struct TP_Modal TP_Modal;
TP_Modal *TP_modal_create(int nvar,int na,int nb,int np,
                          const int *ncols,int *const *cols,double *const *J);
/* Fast vacuum BY construction; U contains the CURRENT physical correction
 * u=(A-1)V at each native node, in Index order. Parameters are process-global. */
TP_Modal *TP_modal_create_analytic(int na,int nb,int np,const double *U);
int TP_modal_solve(TP_Modal *,const double *rhs,double *solution);
void TP_modal_destroy(TP_Modal *);
/* Read-only borrowed factors for the execution adapter, lifetime of m.
 * Block layout [frequency,polar,radial_row,radial_column]. */
int TP_modal_factors(const TP_Modal *, const double **lu,const double **transfer,
                     const double **lower,const double **row_scale,
                     const double **forward,const size_t **permutation);
int TP_modal_shape(const TP_Modal *,int out[3]);
#ifdef __cplusplus
}
#endif
#endif
