/* Benchmark-only wrappers: identical equations and stopping rules, monotonic
 * timing for Newton and its FD setup/preconditioner/spectral JVP phases.
 * Compile instead of TP_Newton.o. */
#include "TwoPunctures.h"
#include "TP_Modal.h"
#include <time.h>
static double elapsed_newton,elapsed_fd_setup,elapsed_modal_setup,elapsed_modal_apply,elapsed_jvp,elapsed_residual;
static double stamp(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static void timed_fd(int,int,int,int,derivs*,int*,int**,double**);
static TP_Modal*timed_modal(int,int,int,const double*);
static int timed_apply(TP_Modal*,const double*,double*);
static void timed_jvp(int,int,int,int,derivs*,double*,derivs*);
static void timed_residual(int,int,int,int,derivs*,double*,derivs*);
#define Newton benchmark_Newton_impl
#define SetMatrix_JFD timed_fd
#define TP_modal_create_analytic timed_modal
#define TP_modal_solve timed_apply
#define J_times_dv timed_jvp
#define F_of_v timed_residual
#include "../src/TP_Newton.c"
#undef Newton
#undef SetMatrix_JFD
#undef TP_modal_create_analytic
#undef TP_modal_solve
#undef J_times_dv
#undef F_of_v
static void timed_fd(int nv,int a,int b,int p,derivs*u,int*nc,int**cols,double**J){double t=stamp();SetMatrix_JFD(nv,a,b,p,u,nc,cols,J);elapsed_fd_setup+=stamp()-t;}
static TP_Modal*timed_modal(int a,int b,int p,const double*u){double t=stamp();TP_Modal*m=TP_modal_create_analytic(a,b,p,u);elapsed_modal_setup+=stamp()-t;return m;}
static int timed_apply(TP_Modal*m,const double*b,double*x){double t=stamp();int status=TP_modal_solve(m,b,x);elapsed_modal_apply+=stamp()-t;return status;}
static void timed_jvp(int nv,int a,int b,int p,derivs*dv,double*out,derivs*u){double t=stamp();J_times_dv(nv,a,b,p,dv,out,u);elapsed_jvp+=stamp()-t;}
static void timed_residual(int nv,int a,int b,int p,derivs*v,double*out,derivs*u){double t=stamp();F_of_v(nv,a,b,p,v,out,u);elapsed_residual+=stamp()-t;}
void Newton(int nv,int a,int b,int p,derivs*v,double tol,int cap){double t=stamp();benchmark_Newton_impl(nv,a,b,p,v,tol,cap);elapsed_newton+=stamp()-t;}
double benchmark_newton_seconds(void){return elapsed_newton;}
void benchmark_phase_seconds(double out[4]){out[0]=elapsed_fd_setup;out[1]=elapsed_modal_setup;out[2]=elapsed_modal_apply;out[3]=elapsed_jvp;}
double benchmark_native_residual_seconds(void){return elapsed_residual;}
