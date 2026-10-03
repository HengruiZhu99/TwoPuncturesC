/* TP_Newton.c */

#include "TwoPunctures.h"
#include "PunctureExecution.h"
#include "TP_LineCache.h"
#include "TP_Modal.h"
#include "PunctureKrylov.h"

static TP_SolverStats solver_stats;
static TP_ExecutionStats execution_stats;
void TP_solver_reset_statistics(void){memset(&solver_stats,0,sizeof(solver_stats));memset(&execution_stats,0,sizeof(execution_stats));}
void TP_solver_get_statistics(TP_SolverStats*out){if(out)*out=solver_stats;}
void TP_solver_get_execution_statistics(TP_ExecutionStats*out){if(out)*out=execution_stats;}

static int
solve_linear (int const nvar, int const n1, int const n2, int const n3,
          derivs* v, derivs* dv,
          int const output, int const itmax, double const tol,
          double * restrict const normres,void *execution_workspace,
          double *frozen_F,derivs *frozen_u);
static double
norm_inf (double const * restrict const F,
          int const ntotal);
static void
relax (double * restrict const dv,
       int const nvar, int const n1, int const n2, int const n3,
       double const * restrict const rhs,
       int * ncols,
       int ** cols,
       double ** JFD, TP_LineCache *cache);
static void
resid (double * restrict const res,
       int const ntotal,
       double const * restrict const dv,
       double const * restrict const rhs,
       int * ncols,
       int ** cols,
       double ** JFD);
static void
LineRelax_al (double * restrict const dv,
              int const j, int const k, int const nvar,
              int const n1, int const n2, int const n3,
	      double const * restrict const rhs,
              int * ncols,
              int ** cols,
              double ** JFD);
static void
LineRelax_be (double * restrict const dv,
              int const i, int const k, int const nvar,
              int const n1, int const n2, int const n3,
	      double const * restrict const rhs,
              int * ncols,
              int ** cols,
              double ** JFD);

/* --------------------------------------------------------------------------*/
static double
norm_inf (double const * restrict const F,
          int const ntotal)
{
  double dmax = -1;
#ifdef TP_OMP
#pragma omp parallel
#endif
  {
    double dmax1 = -1;
#ifdef TP_OMP
#pragma omp for
#endif
    for (int j = 0; j < ntotal; j++)
      if (fabs (F[j]) > dmax1)
        dmax1 = fabs (F[j]);
#ifdef TP_OMP
#pragma omp critical
#endif
    if (dmax1 > dmax)
      dmax = dmax1;
  }
  return dmax;
}

/* --------------------------------------------------------------------------*/
static void
resid (double * restrict const res,
       int const ntotal,
       double const * restrict const dv,
       double const * restrict const rhs,
       int * ncols,
       int ** cols,
       double ** JFD)
{
#ifdef TP_OMP
#pragma omp parallel for
#endif
  for (int i = 0; i < ntotal; i++)
  {
    double JFDdv_i = 0;
    for (int m = 0; m < ncols[i]; m++)
      JFDdv_i += JFD[i][m] * dv[cols[i][m]];
    res[i] = rhs[i] - JFDdv_i;
  }
}

/* -------------------------------------------------------------------------*/
static void
LineRelax_al (double * restrict const dv,
              int const j, int const k, int const nvar,
              int const n1, int const n2, int const n3,
	      double const * restrict const rhs,
              int * ncols,
              int ** cols,
              double ** JFD)
{
  int i, m, Ic, Ip, Im, col, ivar;

  gsl_vector *diag = gsl_vector_alloc(n1);
  gsl_vector *e = gsl_vector_alloc(n1-1); /* above diagonal */
  gsl_vector *f = gsl_vector_alloc(n1-1); /* below diagonal */
  gsl_vector *b = gsl_vector_alloc(n1);   /* rhs */
  gsl_vector *x = gsl_vector_alloc(n1);   /* solution vector */

  for (ivar = 0; ivar < nvar; ivar++)
  {
    gsl_vector_set_zero(diag);
    gsl_vector_set_zero(e);
    gsl_vector_set_zero(f);
    for (i = 0; i < n1; i++)
    {
      Ip = Index (ivar, i + 1, j, k, nvar, n1, n2, n3);
      Ic = Index (ivar, i, j, k, nvar, n1, n2, n3);
      Im = Index (ivar, i - 1, j, k, nvar, n1, n2, n3);
      gsl_vector_set(b,i,rhs[Ic]);
      for (m = 0; m < ncols[Ic]; m++)
      {
	col = cols[Ic][m];
	if (col != Ip && col != Ic && col != Im)
          *gsl_vector_ptr(b, i) -= JFD[Ic][m] * dv[col];
	else
	{
	  if (col == Im && i > 0)
            gsl_vector_set(f,i-1,JFD[Ic][m]);
	  if (col == Ic)
            gsl_vector_set(diag,i,JFD[Ic][m]);
	  if (col == Ip && i < n1-1)
            gsl_vector_set(e,i,JFD[Ic][m]);
	}
      }
    }
    gsl_linalg_solve_tridiag(diag, e, f, b, x);
    for (i = 0; i < n1; i++)
    {
      Ic = Index (ivar, i, j, k, nvar, n1, n2, n3);
      dv[Ic] = gsl_vector_get(x, i);
    }
  }

  gsl_vector_free(diag);
  gsl_vector_free(e);
  gsl_vector_free(f);
  gsl_vector_free(b);
  gsl_vector_free(x);
}

/* --------------------------------------------------------------------------*/
static void
LineRelax_be (double * restrict const dv,
              int const i, int const k, int const nvar,
              int const n1, int const n2, int const n3,
	      double const * restrict const rhs,
              int * ncols,
              int ** cols,
              double ** JFD)
{
  int j, m, Ic, Ip, Im, col, ivar;

  gsl_vector *diag = gsl_vector_alloc(n2);
  gsl_vector *e = gsl_vector_alloc(n2-1); /* above diagonal */
  gsl_vector *f = gsl_vector_alloc(n2-1); /* below diagonal */
  gsl_vector *b = gsl_vector_alloc(n2);   /* rhs */
  gsl_vector *x = gsl_vector_alloc(n2);   /* solution vector */

  for (ivar = 0; ivar < nvar; ivar++)
  {
    gsl_vector_set_zero(diag);
    gsl_vector_set_zero(e);
    gsl_vector_set_zero(f);
    for (j = 0; j < n2; j++)
    {
      Ip = Index (ivar, i, j + 1, k, nvar, n1, n2, n3);
      Ic = Index (ivar, i, j, k, nvar, n1, n2, n3);
      Im = Index (ivar, i, j - 1, k, nvar, n1, n2, n3);
      gsl_vector_set(b,j,rhs[Ic]);
      for (m = 0; m < ncols[Ic]; m++)
      {
	col = cols[Ic][m];
	if (col != Ip && col != Ic && col != Im)
          *gsl_vector_ptr(b, j) -= JFD[Ic][m] * dv[col];
	else
	{
	  if (col == Im && j > 0)
            gsl_vector_set(f,j-1,JFD[Ic][m]);
	  if (col == Ic)
            gsl_vector_set(diag,j,JFD[Ic][m]);
	  if (col == Ip && j < n2-1)
            gsl_vector_set(e,j,JFD[Ic][m]);
	}
      }
    }
    gsl_linalg_solve_tridiag(diag, e, f, b, x);
    for (j = 0; j < n2; j++)
    {
      Ic = Index (ivar, i, j, k, nvar, n1, n2, n3);
      dv[Ic] = gsl_vector_get(x, j);
    }
  }
  gsl_vector_free(diag);
  gsl_vector_free(e);
  gsl_vector_free(f);
  gsl_vector_free(b);
  gsl_vector_free(x);
}

/* --------------------------------------------------------------------------*/
static void
relax (double * restrict const dv,
       int const nvar, int const n1, int const n2, int const n3,
       double const * restrict const rhs,
       int * ncols,
       int ** cols,
       double ** JFD, TP_LineCache *cache)
{
  int i, j, k, n;

  for (k = 0; k < n3; k = k + 2)
  {
    for (n = 0; n < N_PlaneRelax; n++)
    {
#ifdef TP_OMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (i = 2; i < n1; i = i + 2)
	if (cache) TP_line_solve(cache->be + nvar*(i+n1*k), nvar, dv, rhs, ncols, cols, JFD);
        else LineRelax_be (dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#ifdef TP_OMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (i = 1; i < n1; i = i + 2)
	if (cache) TP_line_solve(cache->be + nvar*(i+n1*k), nvar, dv, rhs, ncols, cols, JFD);
        else LineRelax_be (dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#ifdef TP_OMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (j = 1; j < n2; j = j + 2)
	if (cache) TP_line_solve(cache->al + nvar*(j+n2*k), nvar, dv, rhs, ncols, cols, JFD);
        else LineRelax_al (dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#ifdef TP_OMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (j = 0; j < n2; j = j + 2)
	if (cache) TP_line_solve(cache->al + nvar*(j+n2*k), nvar, dv, rhs, ncols, cols, JFD);
        else LineRelax_al (dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
    }
  }
  for (k = 1; k < n3; k = k + 2)
  {
    for (n = 0; n < N_PlaneRelax; n++)
    {
#ifdef TP_OMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (i = 0; i < n1; i = i + 2)
	if (cache) TP_line_solve(cache->be + nvar*(i+n1*k), nvar, dv, rhs, ncols, cols, JFD);
        else LineRelax_be (dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#ifdef TP_OMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (i = 1; i < n1; i = i + 2)
	if (cache) TP_line_solve(cache->be + nvar*(i+n1*k), nvar, dv, rhs, ncols, cols, JFD);
        else LineRelax_be (dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#ifdef TP_OMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (j = 1; j < n2; j = j + 2)
	if (cache) TP_line_solve(cache->al + nvar*(j+n2*k), nvar, dv, rhs, ncols, cols, JFD);
        else LineRelax_al (dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#ifdef TP_OMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (j = 0; j < n2; j = j + 2)
	if (cache) TP_line_solve(cache->al + nvar*(j+n2*k), nvar, dv, rhs, ncols, cols, JFD);
        else LineRelax_al (dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
    }
  }
}

/* --------------------------------------------------------------------------*/
void
TestRelax (int nvar, int n1, int n2, int n3, derivs *v,
	   double *dv)
{
  int ntotal = n1 * n2 * n3 * nvar, **cols, *ncols,
    maxcol = StencilSize * nvar, j;
  double *F, *res, **JFD;
  derivs *u;

  F = dvector (0, ntotal - 1);
  res = dvector (0, ntotal - 1);
  allocate_derivs (&u, ntotal);

  JFD = dmatrix (0, ntotal - 1, 0, maxcol - 1);
  cols = imatrix (0, ntotal - 1, 0, maxcol - 1);
  ncols = ivector (0, ntotal - 1);

  F_of_v (nvar, n1, n2, n3, v, F, u);

  SetMatrix_JFD (nvar, n1, n2, n3, u, ncols, cols, JFD);
  /* Factors belong to this fixed JFD only; no reuse across Newton steps. */
  TP_LineCache *cache = TP_cache_create(nvar, n1, n2, n3, ncols, cols, JFD);

  for (j = 0; j < ntotal; j++)
    dv[j] = 0;
  resid (res, ntotal, dv, F, ncols, cols, JFD);
  printf ("Before: |F|=%20.15e\n", (double) norm1 (res, ntotal));
  fflush(stdout);
  for (j = 0; j < NRELAX; j++)
  {
    relax (dv, nvar, n1, n2, n3, F, ncols, cols, JFD, cache);	/* solves JFD*sh = s*/
    if (j % Step_Relax == 0)
    {
      resid (res, ntotal, dv, F, ncols, cols, JFD);
      printf ("j=%d\t |F|=%20.15e\n", j, (double) norm1 (res, ntotal));
      fflush(stdout);
    }
  }

  resid (res, ntotal, dv, F, ncols, cols, JFD);
  printf ("After: |F|=%20.15e\n", (double) norm1 (res, ntotal));
  fflush(stdout);

  free_dvector (F, 0, ntotal - 1);
  free_dvector (res, 0, ntotal - 1);
  free_derivs (u);

  TP_cache_destroy(cache);
  if(JFD)free_dmatrix (JFD, 0, ntotal - 1, 0, maxcol - 1);
  if(cols)free_imatrix (cols, 0, ntotal - 1, 0, maxcol - 1);
  if(ncols)free_ivector (ncols, 0, ntotal - 1);
}

/* --------------------------------------------------------------------------*/
typedef struct {
 int nvar,n1,n2,n3,total;derivs *u,*direction;
 int *ncols,**cols;double **JFD;TP_LineCache*cache;TP_Modal*modal;
} BYLinearContext;
static int by_linear_action(void*pointer,const double*input,double*output){
 BYLinearContext*c=pointer;double*owned=c->direction->d0;
 c->direction->d0=(double*)input;
 J_times_dv(c->nvar,c->n1,c->n2,c->n3,c->direction,output,c->u);
 c->direction->d0=owned;return 0;
}
static int by_linear_precondition(void*pointer,const double*input,double*output){
 BYLinearContext*c=pointer;memset(output,0,(size_t)c->total*sizeof(double));
 if(c->modal){int status=TP_modal_solve(c->modal,input,output);if(status)solver_stats.modal_failures++;return status;}
 for(int sweep=0;sweep<NRELAX;sweep++){
  relax(output,c->nvar,c->n1,c->n2,c->n3,input,c->ncols,c->cols,c->JFD,c->cache);
  solver_stats.relaxation_sweeps++;
 }
 return 0;
}
static void by_linear_monitor(void*pointer,int iteration,double residual,double alpha,double beta,double omega){
 const char*name=pointer;
 if(iteration==0)printf("%s: %5d  %10.3e\n",name,iteration,residual);
 else printf("%s: %5d  %10.3e  %10.3e  %10.3e  %10.3e\n",name,iteration,residual,alpha,beta,omega);
 fflush(stdout);
}
static double by_dot(const double*a,const double*b,int n){return scalarproduct((double*)a,(double*)b,n);}
static double by_norm(const double*a,int n){return norm2((double*)a,n);}
static int solve_linear(int nvar,int n1,int n2,int n3,derivs*v,derivs*dv,
                         int output,int itmax,double tol,double*normres,void*execution_workspace,
                         double*frozen_F,derivs*frozen_u){
 BYLinearContext c={0};c.nvar=nvar;c.n1=n1;c.n2=n2;c.n3=n3;c.total=nvar*n1*n2*n3;
 const int modal_requested=params_get_int("TP_preconditioner"),method=params_get_int("TP_krylov_solver");
 const int strict=modal_requested||params_get_int("TP_linear_relative")||method!=PK_BICGSTAB;
 double*F;
 if(execution_workspace){F=frozen_F;c.u=frozen_u;}
 else{F=dvector(0,c.total-1);allocate_derivs(&c.u,c.total);allocate_derivs(&c.direction,c.total);
  F_of_v(nvar,n1,n2,n3,v,F,c.u);}
 int failed=0;
 if(modal_requested){
  c.modal=nvar==1?TP_modal_create_analytic(n1,n2,n3,c.u->d0):NULL;
  if(c.modal)solver_stats.modal_factorizations++;else{solver_stats.modal_failures++;failed=1;}
 }else{
  c.JFD=dmatrix(0,c.total-1,0,StencilSize*nvar-1);c.cols=imatrix(0,c.total-1,0,StencilSize*nvar-1);c.ncols=ivector(0,c.total-1);
  SetMatrix_JFD(nvar,n1,n2,n3,c.u,c.ncols,c.cols,c.JFD);
  c.cache=TP_cache_create(nvar,n1,n2,n3,c.ncols,c.cols,c.JFD);
 }
 const char*name=method==PK_GMRES?"gmres":"bicgstab";
 if(output){printf("%s:  itmax %d, tol %e\n",name,itmax,tol);fflush(stdout);}
 PK_Options options={method,itmax,params_get_int("TP_krylov_restart"),strict,method==PK_BICGSTAB,0,tol,method==PK_BICGSTAB?by_dot:NULL,method==PK_BICGSTAB?by_norm:NULL};
 PK_Result result={0};
 result.true_residual=result.relative_residual=NAN;
 int status;
#ifdef PUNCTURES_KOKKOS
 if(params_get_int("TP_execution_backend"))status=failed?PK_CALLBACK:Puncture_kokkos_BY_linear(execution_workspace,c.u->d0,c.modal,F,dv->d0,&options,&result,output?by_linear_monitor:NULL,(void*)name);
 else
#endif
 status=failed?PK_CALLBACK:PK_solve(c.total,F,dv->d0,&options,by_linear_action,by_linear_precondition,&c,output?by_linear_monitor:NULL,(void*)name,&result);
 *normres=result.recurrence_residual;
 solver_stats.krylov_iterations+=result.iterations;solver_stats.jvp_applications+=result.operator_calls;
 solver_stats.preconditioner_applications+=result.preconditioner_calls;solver_stats.last_linear_target=tol;
 if(strict){
  solver_stats.last_true_linear_residual=result.true_residual;solver_stats.last_relative_linear_residual=result.relative_residual;
  if(status!=PK_SUCCESS)failed=1;
  if(output)printf("linear_true: %.17e relative %.17e target %.17e\n",result.true_residual,result.relative_residual,tol);
  if(status!=PK_SUCCESS&&output)printf("linear_failure: %d %s\n",status,PK_status_string(status));
 }
 if(failed)solver_stats.linear_failures++;
 TP_modal_destroy(c.modal);TP_cache_destroy(c.cache);
 if(c.JFD)free_dmatrix(c.JFD,0,c.total-1,0,StencilSize*nvar-1);
 if(c.cols)free_imatrix(c.cols,0,c.total-1,0,StencilSize*nvar-1);
 if(c.ncols)free_ivector(c.ncols,0,c.total-1);
 if(!execution_workspace){free_dvector(F,0,c.total-1);free_derivs(c.u);free_derivs(c.direction);}
 // Legacy Newton ignores a recurrence failure, as before; new modes fail closed.
 return failed?-(100+status):status==PK_SUCCESS?0:result.iterations;
}

static int newton_residual(int nvar,int n1,int n2,int n3,derivs*v,double*F,derivs*u,void*workspace,double tol,int*reference_polishing){
#ifdef PUNCTURES_KOKKOS
 if(workspace&&!*reference_polishing){
  int status=Puncture_kokkos_BY_residual(workspace,v,F,u);if(status)return status;
  execution_stats.device_residual_calls++;execution_stats.device_candidate_linf=norm_inf(F,nvar*n1*n2*n3);
  if(execution_stats.device_candidate_linf>tol)return 0;
  // A cached matrix and a transform can differ at roundoff near the map
  // foci. Confirm the original residual before accepting outer convergence,
  // then retain it for any remaining refinement of this same iterate.
  *reference_polishing=1;
  execution_stats.original_confirmation_calls++;
 }
#endif
 F_of_v(nvar,n1,n2,n3,v,F,u);
 if(workspace){
  const double*arrays[21]={F,v->d0,v->d1,v->d2,v->d3,v->d11,v->d12,v->d13,v->d22,v->d23,v->d33,u->d0,u->d1,u->d2,u->d3,u->d11,u->d12,u->d13,u->d22,u->d23,u->d33};
  for(int a=0;a<21;a++)for(int p=0;p<nvar*n1*n2*n3;p++)if(!isfinite(arrays[a][p]))return PK_NONFINITE;
  execution_stats.original_confirmation_linf=norm_inf(F,nvar*n1*n2*n3);
 }return 0;
}

/* -------------------------------------------------------------------*/
void
Newton (int const nvar, int const n1, int const n2, int const n3,
	derivs *v,
        double const tol, int const itmax)
{
  int verbose = params_get_int("verbose");
  
  int ntotal = n1 * n2 * n3 * nvar, ii, it;
  double *F, dmax, normres;
  derivs *u, *dv;
  
  F = dvector (0, ntotal - 1);
  allocate_derivs (&dv, ntotal);
  allocate_derivs (&u, ntotal);
  void *execution_workspace=NULL;
  int reference_polishing=0;
#ifdef PUNCTURES_KOKKOS
  if(params_get_int("TP_execution_backend")){
    if(nvar==1&&!params_get_int("do_residuum_debug_output"))execution_workspace=Puncture_kokkos_BY_create(n1,n2,n3);
    if(!execution_workspace){solver_stats.linear_failures++;free_dvector(F,0,ntotal-1);free_derivs(dv);free_derivs(u);return;}
    execution_stats.workspace_creations++;
  }
#endif
  
  /*         TestRelax(nvar, n1, n2, n3, v, dv->d0); */
  it = 0;
  dmax = 1;
  while (dmax > tol && it < itmax)
    {
      if (it == 0)
	{
	  if(newton_residual(nvar,n1,n2,n3,v,F,u,execution_workspace,tol,&reference_polishing)){solver_stats.linear_failures++;break;}
	  dmax = norm_inf (F, ntotal);
	}
#ifdef TP_OMP
#pragma omp parallel for
#endif
      for (int j = 0; j < ntotal; j++)
	dv->d0[j] = 0;
      
      if(verbose){
	printf ("Newton: it=%d \t |F|=%e\n", it, (double)dmax);
	printf ("bare mass: mp=%g \t mm=%g\n",
		params_get_real("par_m_plus"), params_get_real("par_m_minus"));
      }
      
      fflush(stdout);
      double linear_target=dmax*params_get_real("TP_linear_rtol");
      if(params_get_int("TP_linear_relative"))linear_target=norm2(F,ntotal)*params_get_real("TP_linear_rtol");
      ii = solve_linear (nvar, n1, n2, n3, v, dv, verbose, params_get_int("TP_krylov_maxit"), linear_target, &normres,execution_workspace,F,u);
      if((params_get_int("TP_preconditioner")||params_get_int("TP_linear_relative")||params_get_int("TP_krylov_solver")!=PK_BICGSTAB)&&ii<0){
        if(verbose)printf("Newton linear solve failed: %d\n",ii);break;
      }
      solver_stats.newton_iterations++;
      if(execution_workspace&&reference_polishing)execution_stats.polishing_steps++;
#ifdef TP_OMP
#pragma omp parallel for
#endif
      for (int j = 0; j < ntotal; j++)
	v->d0[j] -= dv->d0[j];
      if(newton_residual(nvar,n1,n2,n3,v,F,u,execution_workspace,tol,&reference_polishing)){solver_stats.linear_failures++;break;}
      dmax = norm_inf (F, ntotal);
      it += 1;
    }
  if (itmax==0)
    {
      if(newton_residual(nvar,n1,n2,n3,v,F,u,execution_workspace,tol,&reference_polishing))solver_stats.linear_failures++;
      dmax = norm_inf (F, ntotal);
    }
  
  if(verbose)
    printf ("Newton: it=%d \t |F|=%e \n", it, (double)dmax);
  
  fflush(stdout);
  
  free_dvector (F, 0, ntotal - 1);
  free_derivs (dv);
  free_derivs (u);
#ifdef PUNCTURES_KOKKOS
  Puncture_kokkos_BY_destroy(execution_workspace);
#endif
}
