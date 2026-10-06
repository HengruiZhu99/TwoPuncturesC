#ifndef HISPID_H
#define HISPID_H
#include "PunctureKrylov.h"
#include "PunctureExecution.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Version 1. Units G=c=1, signature -+++, K_ij=-1/2 L_n gamma_ij.
 * Arrays are row-major; symmetric tensors use FULL 3x3 storage.
 * Seeds: rest Kerr mass, rest Cartesian spin S, lab velocity v (not ADM P).
 * |S|<mass^2, |v|<1. mass=0 disables a hole. */
typedef struct {
  double mass, center[3], spin[3], velocity[3];
} HiSpID_Hole;
/* Independent extents; aggregate host/device memory guards still apply. */
enum {HISPID_MAX_RADIAL_POINTS=512,HISPID_MAX_POLAR_POINTS=512,
      HISPID_MAX_AZIMUTHAL_POINTS=256};
typedef struct {
  HiSpID_Hole hole[2];
  int n[3];
  int conformal_choice; /* 0: QI scalar, 1: det(gamma)^(1/12) */
  int inner_flatten;    /* Eq27/28 modified correction operators */
  double omega[2];      /* cross attenuation widths; <=0 disables f */
  int attenuation_power;
  double inner_min[2], inner_max[2]; /* <=0 disables g for that hole */
  double far_radius;    /* <=0 disables F */
  double tolerance;
  int max_newton, max_krylov, krylov_restart;
 int memory_limit_mib; /* conservative aggregate budget; default2048, max65536 */
} HiSpID_Config;
typedef struct {
  double gamma[9], Kij[9], psi, conformal_metric[9], Atilde[9];
  double mean_curvature, correction[4], attenuation;
} HiSpID_Point;
typedef struct {
  int converged, newton_iterations, krylov_iterations, npoints;
  double scaled_linf[4], unscaled_linf[4], seconds;
} HiSpID_Diagnostics;
typedef struct HiSpID_Data HiSpID_Data;

/* ABI-safe family selection: existing configuration layouts/entry points
 * retain their QI meaning. The family applies to both active holes. */
enum HiSpID_SeedFamily {HISPID_SEED_QI=0,HISPID_SEED_TRUMPET_R0_M=1};
HiSpID_Data *HiSpID_create_with_seed_family(const HiSpID_Config *,int family,
                                          int execution,int geometry_execution,
                                          int sampler_only);
int HiSpID_seed_family(const HiSpID_Data *);
int HiSpID_seed_with_family(const HiSpID_Hole *,int conformal_choice,int count,
                            const double *xyz,HiSpID_Point *,int execution,int family);
int HiSpID_operators_with_seed_family(const HiSpID_Config *,const double xyz[3],
                                     const double jets[40],double out[5],int family);

void HiSpID_default_config(HiSpID_Config *);
/* Identifies the continuous basis represented by the saved nodal unknowns. */
const char *HiSpID_unknown_parameterization(void);
/* Fixed build maps [radial_stretch,angular_stretch], also encoded in the
 * exact basis identifier for nondefault builds. Return0 on success. */
int HiSpID_collocation_maps(double out[2]);
/* Names the build's positive residual/JVP/preconditioner row scaling. */
const char *HiSpID_residual_scaling(void);
/* No global parameter changes. Contexts own their configuration and caches.
 * Both map centers must be distinct, even when one mass is zero.
 * Return NULL on invalid configuration/allocation/geometry failure. */
HiSpID_Data *HiSpID_create(const HiSpID_Config *);
/* Explicit execution choice; reference is unchanged. Kokkos uses this
 * image's execution space and parallel host geometry setup when available. */
HiSpID_Data *HiSpID_create_with_execution(const HiSpID_Config *,int execution);
/* Explicit setup choice: geometry_execution=0 retains the long-double host
 * reference. Value1 builds seed/coordinate/derivative/operator caches in
 * the selected Kokkos execution space, using double precision. It requires
 * execution=PUNCTURE_KOKKOS. Existing constructors retain geometry_execution=0. */
HiSpID_Data *HiSpID_create_with_geometry(const HiSpID_Config *,int execution,
                                       int geometry_execution);
/* Setup phase wall times; coefficient_seconds accumulates lazy transforms.
 * scalar_digits is the geometry scalar's binary precision. Initialize
 * struct_size=sizeof(HiSpID_SetupStatistics) before querying. */
typedef struct {
 int struct_size,geometry_execution,scalar_digits;
 double spectral_seconds,geometry_seconds,coefficient_seconds;
} HiSpID_SetupStatistics;
int HiSpID_setup_statistics(const HiSpID_Data *,HiSpID_SetupStatistics *);
/* Sampling-only context for loading saved unknowns: no collocation geometry,
 * derivative workspace, or Newton/Krylov allocation. It supports sampling,
 * charges and get/set_unknowns, and rejects solve/residual/JVP operations. */
HiSpID_Data *HiSpID_create_sampler(const HiSpID_Config *);
/* 0 converged, 1 iteration/line-search/Krylov failure, -1 invalid context.
 * Failed solves remain sampleable for diagnosis. */
int HiSpID_solve(HiSpID_Data *);
/* Opt-in no-swirl axial sector for exactly coaxial, nonspinning data.
 * Defaults OFF. Only the linear search is projected: Newton convergence and
 * line search still use the full, unprojected equations. Does not change the
 * stored field representation or checkpoint format. Rejects invalid geometry. */
int HiSpID_set_axisymmetric(HiSpID_Data *, int enabled);
/* Optional fixed RHS-relative L2 forcing, 0<rtol<1. The default entry point
 * retains its adaptive forcing. Residual row scaling is fixed by the build. */
int HiSpID_solve_with_forcing(HiSpID_Data *, double rtol);
/* ABI-safe options separate the linear backend from the equation-system config.
 * krylov: PK_GMRES=0, PK_BICGSTAB=1 or PK_LGMRES=2 (PunctureKrylov.h).
 * rtol=0 retains adaptive forcing; positive rtol selects fixed relative L2. */
typedef struct { int struct_size,krylov; double linear_rtol; } HiSpID_SolveOptions;
void HiSpID_default_solve_options(HiSpID_SolveOptions *);
int HiSpID_solve_with_options(HiSpID_Data *, const HiSpID_SolveOptions *);
/* Last options actually accepted by the native solve; -1 before a solve. */
int HiSpID_resolved_solve_options(const HiSpID_Data *,HiSpID_SolveOptions *);
/* Work counters [spectral JVP calls, preconditioner applies]. */
int HiSpID_work_statistics(const HiSpID_Data *, int out[2]);
/* Solves: rows [Newton index, requested relative L2 target,
 * true relative linear residual, Krylov count]. capacity=0/out=NULL queries
 * count; return count, or -1 on invalid buffer/context. */
int HiSpID_linear_history(const HiSpID_Data *, int capacity, double *out);
int HiSpID_diagnostics(const HiSpID_Data *, HiSpID_Diagnostics *);
int HiSpID_sample(HiSpID_Data *, int count, const double *xyz, HiSpID_Point *);
/* Physical lab-frame metric gradients, dgamma[27*p+9*d+3*i+j]=d_d gamma_ij.
 * Correction values and gradients use the same Cartesian C2 modal basis,
 * including its analytic limits at the prolate axes (map foci excluded). */
int HiSpID_sample_with_derivatives(HiSpID_Data *, int count, const double *xyz,
                                   HiSpID_Point *, double *dgamma);
int HiSpID_seed(const HiSpID_Hole *, int conformal_choice,
                int count, const double *xyz, HiSpID_Point *);
/* Optional batch seed evaluator. PUNCTURE_REFERENCE calls the unchanged
 * long-double host API; PUNCTURE_KOKKOS evaluates stable double geometry in
 * this image's execution space. xyz/output are host buffers. Output is valid
 * only on return0: failed device status/finiteness never falls back to host.
 * No solve, attenuation, corrections or derivative jets are exported. */
int HiSpID_seed_with_execution(const HiSpID_Hole *, int conformal_choice,
                              int count, const double *xyz, HiSpID_Point *,
                              int execution);
/* Finite-radius PHYSICAL ADM surface integrals E,P[3],J[3] in lab coordinates,
 * angular momentum about sphere_center. Extrapolation is caller's job. */
int HiSpID_charges(HiSpID_Data *, const double sphere_center[3], double radius,
                  int ntheta, int nphi, double EPJ[7]);
/* Native manufactured-operator hook: input jets [4,10] = value, dx,dy,dz,
 * dxx,dxy,dxz,dyy,dyz,dzz. Outputs scalar Laplacian, vector longitudinal
 * divergence, and scalar curvature, with g=1 and actual background metric. */
int HiSpID_operators(const HiSpID_Config *, const double xyz[3],
                    const double jets[40], double out[5]);
int HiSpID_get_unknowns(const HiSpID_Data *, double *values, int length);
int HiSpID_set_unknowns(HiSpID_Data *, const double *values, int length);
int HiSpID_residual(HiSpID_Data *, const double *values, double *residual);
int HiSpID_jvp(HiSpID_Data *, const double *base, const double *direction,
               double *out);
/* npoints collocation samples: xyz[3*npoints], g[npoints], psi[npoints],
 * HM[4*npoints] = [-8 FH/psi^5, FM^x/psi^10,...]. These equal physical
 * constraints ONLY where g=1. No independent physical verification implied. */
int HiSpID_equation_samples(HiSpID_Data *, double *xyz, double *g,
                            double *psi, double *HM);
void HiSpID_destroy(HiSpID_Data *);
const char *HiSpID_last_error(void);

#ifdef __cplusplus
}
#endif
#endif
