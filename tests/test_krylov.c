/* Independent manufactured nonsymmetric matrix controls for both backends. */
#include "PunctureKrylov.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#define N 9
static int checks;
#define CHECK(c) do {checks++;if(!(c)){fprintf(stderr,"failed line %d: %s\n",__LINE__,#c);return 1;}}while(0)
typedef struct { double a[N][N]; int fail,nan,diagonal,pre_fail,pre_nan; } Matrix;
static int apply(void *v,const double *x,double *y){
  Matrix*m=v;if(m->fail)return -1;
  for(int i=0;i<N;i++){y[i]=0;for(int j=0;j<N;j++)y[i]+=m->a[i][j]*x[j];}
  if(m->nan)y[0]=NAN;return 0;
}
static int precondition(void*v,const double*x,double*y){
  Matrix*m=v;if(m->pre_fail)return -1;for(int i=0;i<N;i++)y[i]=m->diagonal?x[i]/m->a[i][i]:x[i];if(m->pre_nan)y[0]=NAN;return 0;
}
static double residual(Matrix*m,const double*b,const double*x){
  double sum=0;for(int i=0;i<N;i++){double r=b[i];for(int j=0;j<N;j++)r-=m->a[i][j]*x[j];sum+=r*r;}return sqrt(sum);
}
static int null_operator(void*c,const double*x,double*y){(void)c;y[0]=0;y[1]=x[1];return 0;}
static int huge_precondition(void*c,const double*x,double*y){(void)c;y[0]=1e308*x[1];y[1]=x[1];return 0;}
int main(void){
  Matrix m={0};double exact[N],rhs[N],x[N];PK_Result result;
  for(int i=0;i<N;i++){exact[i]=sin(i+.25);for(int j=0;j<N;j++)m.a[i][j]=i==j?4+i*.1:(j==i+1?-.8:(i==j+2?.35:0));}
  apply(&m,exact,rhs);
  for(int method=0;method<2;method++)for(int diagonal=0;diagonal<2;diagonal++)for(int guess=0;guess<2;guess++){
    PK_Options o={method,200,3,1,0,0,1e-11,NULL,NULL};m.diagonal=diagonal;
    for(int i=0;i<N;i++)x[i]=guess?.3*cos(i):0;
    CHECK(PK_solve(N,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_SUCCESS);
    CHECK(residual(&m,rhs,x)<=o.absolute_tolerance);
    CHECK(result.true_residual<=o.absolute_tolerance);
    CHECK(result.operator_calls>result.iterations&&result.preconditioner_calls>=result.iterations);
    for(int i=0;i<N;i++)CHECK(fabs(x[i]-exact[i])<1e-10);
  }
  for(int method=0;method<2;method++){
    memset(&m,0,sizeof(m));for(int i=0;i<N;i++)for(int j=0;j<N;j++)m.a[i][j]=i==j?4+i*.1:(j==i+1?-.8:(i==j+2?.35:0));
    PK_Options o={method,1,4,1,0,0,1e-15,NULL,NULL};memset(x,0,sizeof(x));m.diagonal=0;
    CHECK(PK_solve(N,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_LIMIT);
    CHECK(result.iterations==1&&residual(&m,rhs,x)>o.absolute_tolerance);
    m.pre_fail=1;memset(x,0,sizeof(x));CHECK(PK_solve(N,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_CALLBACK);m.pre_fail=0;
    m.pre_nan=1;CHECK(PK_solve(N,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_NONFINITE);m.pre_nan=0;
    m.fail=1;CHECK(PK_solve(N,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_CALLBACK);m.fail=0;
    m.nan=1;CHECK(PK_solve(N,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_NONFINITE);m.nan=0;
    memset(&m,0,sizeof(m));memset(x,0,sizeof(x));
    CHECK(PK_solve(N,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_BREAKDOWN);
    for(int i=0;i<N;i++)m.a[i][i]=1;
    o.max_iterations=20;o.absolute_tolerance=1e-12;
    CHECK(PK_solve(N,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_SUCCESS);
    CHECK(result.iterations==1&&residual(&m,rhs,x)<=o.absolute_tolerance);
    double zero[N]={0};memset(x,0,sizeof(x));
    CHECK(PK_solve(N,zero,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_SUCCESS);
    CHECK(result.iterations==0&&result.preconditioner_calls==0&&result.true_residual==0);
    CHECK(result.operator_calls==(method==PK_GMRES?0:2));
    o.method=2;CHECK(PK_solve(N,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_INVALID);o.method=method;
    o.restart=4097;CHECK(PK_solve(N,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_INVALID);o.restart=4;
    o.verify_true_residual=2;CHECK(PK_solve(N,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_INVALID);o.verify_true_residual=1;
    CHECK(PK_solve(N,rhs,rhs,&o,apply,precondition,&m,NULL,NULL,&result)==PK_INVALID);
    CHECK(PK_solve(0,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_INVALID);
    o.absolute_tolerance=NAN;CHECK(PK_solve(N,rhs,x,&o,apply,precondition,&m,NULL,NULL,&result)==PK_INVALID);
    CHECK(PK_solve(N,rhs,x,&o,NULL,precondition,&m,NULL,NULL,&result)==PK_INVALID);
  }
  {double b[2]={0,2},z[2]={0};PK_Options o={PK_GMRES,10,2,1,0,0,1e-10,NULL,NULL};
   CHECK(PK_solve(2,b,z,&o,null_operator,huge_precondition,NULL,NULL,NULL,&result)==PK_NONFINITE);}
  printf("shared Krylov independent controls: %d passed\n",checks);return 0;
}
