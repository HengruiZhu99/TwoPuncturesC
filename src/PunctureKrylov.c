#include "PunctureKrylov.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
static double dot(const double*a,const double*b,int n){double q=0;for(int i=0;i<n;i++)q+=a[i]*b[i];return q;}
static double length(const double*a,int n){return sqrt(dot(a,a,n));}
static double inner(const PK_Options*o,const double*a,const double*b,int n){return o->dot?o->dot(a,b,n):dot(a,b,n);}
static double norm(const PK_Options*o,const double*a,int n){return o->norm?o->norm(a,n):length(a,n);}
static int finite_vector(const double*a,int n){for(int i=0;i<n;i++)if(!isfinite(a[i]))return 0;return 1;}
static int action(PK_Apply f,void*c,const double*b,double*x,int n,int*count){
  (*count)++;if(f(c,b,x))return PK_CALLBACK;return finite_vector(x,n)?PK_SUCCESS:PK_NONFINITE;
}
static void *bank(size_t rows,size_t cols,size_t width){
  if(cols&&rows>SIZE_MAX/cols)return NULL;size_t n=rows*cols;
  if(width&&n>SIZE_MAX/width)return NULL;return calloc(n,width);
}
static double *vector_bank(size_t n){return bank(n,1,sizeof(double));}
static void vector_free(double*v){free(v);}
static void vector_zero(double*v,int n){memset(v,0,(size_t)n*sizeof(double));}
static void vector_difference(double*r,const double*b,const double*Ax,int n){for(int i=0;i<n;i++)r[i]=b[i]-Ax[i];}
static void vector_divide(double*out,const double*in,double beta,int n){for(int i=0;i<n;i++)out[i]=in[i]/beta;}
static void vector_subtract(double*w,const double*v,double a,int n){for(int q=0;q<n;q++)w[q]-=a*v[q];}
static void vector_add(double*x,const double*v,double a,int n){for(int q=0;q<n;q++)x[q]+=v[q]*a;}
static void vector_add_left(double*x,const double*v,double a,int n){for(int j=0;j<n;j++)x[j]+=a*v[j];}
static void vector_shadow(double*rt,double*r,const double*b,int n){for(int j=0;j<n;j++)rt[j]=r[j]=b[j]-r[j];}
static void vector_copy(double*p,const double*r,int n){for(int j=0;j<n;j++)p[j]=r[j];}
static void vector_search(double*p,const double*r,const double*v,double beta,double omega,int n){for(int j=0;j<n;j++)p[j]=r[j]+beta*(p[j]-omega*v[j]);}
static void vector_step(double*s,const double*r,const double*v,double alpha,int n){for(int j=0;j<n;j++)s[j]=r[j]-alpha*v[j];}
static void vector_bicg_update(double*x,double*r,const double*s,const double*t,const double*ph,const double*sh,double alpha,double omega,int n){for(int j=0;j<n;j++){x[j]+=alpha*ph[j]+omega*sh[j];r[j]=s[j]-omega*t[j];}}
#include "PunctureKrylov_controller.inc"
int PK_solve(int n,const double*b,double*x,const PK_Options*o,PK_Apply A,PK_Apply M,void*c,PK_Monitor monitor,void*mc,PK_Result*out){return controller_solve(n,b,x,o,A,M,c,monitor,mc,out);}
const char *PK_status_string(int status){
  switch(status){
    case PK_SUCCESS:return "converged";
    case PK_LIMIT:return "iteration limit or true residual above target";
    case PK_BREAKDOWN:return "Krylov breakdown";
    case PK_CALLBACK:return "operator/preconditioner callback failed";
    case PK_INVALID:return "invalid Krylov options or vectors";
    case PK_ALLOCATION:return "Krylov workspace allocation failed";
    case PK_NONFINITE:return "nonfinite Krylov result";
    default:return "unknown Krylov status";
  }
}
