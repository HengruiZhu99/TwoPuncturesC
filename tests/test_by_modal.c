#include "TwoPunctures.h"
#include "TP_Modal.h"
#include <limits.h>
/* Independent physical-grid projection: average each actual JFD stencil over
 * azimuth, translate it around the full circle, and symmetrize opposite shifts.
 * No modal assembly/factor code is used to build this reference matrix. */
static double *project(int na,int nb,int np,int*nc,int**col,double**J){
 int stride=na*nb,N=stride*np;double*A=calloc((size_t)N*N,sizeof(double));if(!A)return NULL;
 for(int k=0;k<np;k++)for(int p=0;p<stride;p++){
  int row=p+stride*k;
  for(int q=0;q<nc[row];q++){
   int c=col[row][q]%stride,shift=col[row][q]/stride-k;double value=J[row][q]/(2*np);
   for(int t=0;t<np;t++){
    int plus=(t+shift+np)%np,minus=(t-shift+np)%np;
    A[(size_t)(p+stride*t)*N+c+stride*plus]+=value;
    A[(size_t)(p+stride*t)*N+c+stride*minus]+=value;
   }
  }
 }
 return A;
}
int main(void){
 TwoPunctures_params_set_default();params_set_int("verbose",0);params_set_real("par_b",3);
 params_set_real("par_m_plus",.6);params_set_real("par_m_minus",.4);
 params_set_real("par_P_plus2",.036);params_set_real("par_S_plus3",.108);params_set_real("par_S_minus1",-.032);
 int checks=0;double largest_solution=0,largest_residual=0;
 for(int fixture=0;fixture<2;fixture++){
  int na=fixture?7:6,nb=fixture?9:8,np=fixture?8:6,N=na*nb*np;
  derivs*v,*u;allocate_derivs(&v,N);allocate_derivs(&u,N);double*F=dvector(0,N-1);
  int*nc=ivector(0,N-1),**cols=imatrix(0,N-1,0,StencilSize-1);double**J=dmatrix(0,N-1,0,StencilSize-1);
  double*known=dvector(0,N-1),*rhs=dvector(0,N-1),*x=dvector(0,N-1);
  for(int refresh=0;refresh<2;refresh++){
   params_set_real("par_m_minus",.4-.02*refresh);
   for(int q=0;q<N;q++)v->d0[q]=1e-4*(refresh+1)*cos(.19*q);
   F_of_v(1,na,nb,np,v,F,u);SetMatrix_JFD(1,na,nb,np,u,nc,cols,J);
   TP_Modal*m=TP_modal_create(1,na,nb,np,nc,cols,J);if(!m){fprintf(stderr,"factor failed\n");return 1;}
   TP_Modal*fast=TP_modal_create_analytic(na,nb,np,u->d0);if(!fast)return 8;
   double*A=project(na,nb,np,nc,cols,J);if(!A)return 2;
   for(int mode=0;mode<np;mode++){
    for(int k=0;k<np;k++)for(int p=0;p<na*nb;p++)known[p+na*nb*k]=(1+.1*cos(.13*p))*cos(2*Pi*mode*k/np+.31);
    for(int row=0;row<N;row++){rhs[row]=0;for(int col=0;col<N;col++)rhs[row]+=A[(size_t)row*N+col]*known[col];}
    if(TP_modal_solve(m,rhs,x))return 3;
    double difference=0,residual=0,scale=0;
    for(int row=0;row<N;row++){
     difference=fmax(difference,fabs(x[row]-known[row]));double value=0;
     for(int col=0;col<N;col++)value+=A[(size_t)row*N+col]*x[col];
     residual=fmax(residual,fabs(value-rhs[row]));scale=fmax(scale,fabs(rhs[row]));
    }
    largest_solution=fmax(largest_solution,difference);largest_residual=fmax(largest_residual,residual/fmax(1,scale));
    if(difference>2e-9||residual/fmax(1,scale)>2e-12){fprintf(stderr,"inverse mismatch %.17e %.17e\n",difference,residual/fmax(1,scale));return 4;}checks++;
    if(TP_modal_solve(fast,rhs,x))return 9;
    difference=0;residual=0;
    for(int row=0;row<N;row++){
     difference=fmax(difference,fabs(x[row]-known[row]));double value=0;
     for(int col=0;col<N;col++)value+=A[(size_t)row*N+col]*x[col];
     residual=fmax(residual,fabs(value-rhs[row]));
    }
    largest_solution=fmax(largest_solution,difference);largest_residual=fmax(largest_residual,residual/fmax(1,scale));
    if(difference>2e-9||residual/fmax(1,scale)>2e-12){fprintf(stderr,"analytic/original mismatch %.17e %.17e\n",difference,residual/fmax(1,scale));return 10;}checks++;
   }
   double saved_rhs=rhs[0];rhs[0]=NAN;
   if(!TP_modal_solve(m,rhs,x)||!TP_modal_solve(fast,rhs,x))return 11;checks+=2;
   rhs[0]=saved_rhs;if(TP_modal_solve(m,rhs,x)||TP_modal_solve(fast,rhs,x))return 12;checks+=2;
   free(A);TP_modal_destroy(m);TP_modal_destroy(fast);
   if(TP_modal_create(2,na,nb,np,nc,cols,J)||TP_modal_create(1,na,nb,np-1,nc,cols,J))return 5;checks+=2;
   double saved=J[0][0];J[0][0]=NAN;if(TP_modal_create(1,na,nb,np,nc,cols,J))return 6;J[0][0]=saved;checks++;
   int*saved_col=cols[0];cols[0]=NULL;if(TP_modal_create(1,na,nb,np,nc,cols,J))return 13;cols[0]=saved_col;checks++;
   double*saved_row=J[0];J[0]=NULL;if(TP_modal_create(1,na,nb,np,nc,cols,J))return 14;J[0]=saved_row;checks++;
   if(TP_modal_create(1,INT_MAX,INT_MAX,INT_MAX-1,nc,cols,J)||TP_modal_create_analytic(INT_MAX,INT_MAX,INT_MAX-1,u->d0))return 15;checks+=2;
   saved=u->d0[0];u->d0[0]=NAN;if(TP_modal_create_analytic(na,nb,np,u->d0))return 16;u->d0[0]=saved;checks++;
   saved=params_get_real("par_b");params_set_real("par_b",1e-300);
   if(TP_modal_create_analytic(na,nb,np,u->d0))return 24;params_set_real("par_b",saved);checks++;
   for(int row=0;row<N;row++){nc[row]=1;cols[row][0]=na*((row/na)%nb)+na*nb*(row/(na*nb));J[row][0]=1;}
   if(TP_modal_create(1,na,nb,np,nc,cols,J))return 17;checks++;
   for(int row=0;row<N;row++)for(int q=0;q<nc[row];q++)J[row][q]=0;
   if(TP_modal_create(1,na,nb,np,nc,cols,J))return 7;checks++;
  }
  free_derivs(v);free_derivs(u);free_dvector(F,0,N-1);free_ivector(nc,0,N-1);free_imatrix(cols,0,N-1,0,StencilSize-1);free_dmatrix(J,0,N-1,0,StencilSize-1);
  free_dvector(known,0,N-1);free_dvector(rhs,0,N-1);free_dvector(x,0,N-1);
 }
 TwoPunctures_params_set_default();
 params_set_int("npoints_A",6);params_set_int("npoints_B",8);params_set_int("npoints_phi",6);
 for(int selector=0;selector<2;selector++){
  const char*name=selector?"TP_linear_relative":"TP_preconditioner";
  params_set_int(name,2);if(TwoPunctures_make_initial_data())return 18;params_set_int(name,0);checks++;
 }
 double invalid[]={NAN,0,1,-.01,INFINITY};
 for(size_t q=0;q<sizeof(invalid)/sizeof(invalid[0]);q++){
  params_set_real("TP_linear_rtol",invalid[q]);if(TwoPunctures_make_initial_data())return 19;checks++;
 }
 params_set_int("TP_preconditioner",1);params_set_real("TP_linear_rtol",1e-30);
 params_set_real("par_b",3);params_set_real("par_m_plus",.6);params_set_real("par_m_minus",.4);
 params_set_real("par_S_plus3",.108);params_set_int("give_bare_mass",0);
 ini_data*failed=TwoPunctures_make_initial_data();if(!failed)return 20;
 TP_SolverStats stats;TP_solver_get_statistics(&stats);
 if(stats.linear_failures!=1||stats.newton_iterations!=0||params_get_real("par_m_plus")!=.6||params_get_real("par_m_minus")!=.4)return 21;
 for(int q=0;q<failed->ntotal;q++)if(failed->v->d0[q]!=0)return 22;
 double res,energy,masses[2];if(TwoPunctures_diagnostics(failed,&res,&energy,masses)!=1)return 23;
 TwoPunctures_finalise(failed);checks+=4;
 printf("BY modal independent cyclic projection/inverse: %d checks, solution max %.17e, residual max %.17e\n",checks,largest_solution,largest_residual);return 0;
}
