/* Actual BY spectral Jacobian and full-lower modal inverse. The inverse
 * oracle is independently projected from the original physical JFD rows. */
#include "PunctureKokkos.hpp"
#include "TP_Modal.h"
extern "C" {
#include "TwoPunctures.h"
}
#include <cstdio>
using namespace puncture;
static int checks;
static double forbidden_dot(const double*,const double*,int){std::abort();}
static double forbidden_norm(const double*,int){std::abort();}
#define CHECK(c) do{checks++;if(!(c)){fprintf(stderr,"line%d: %s\n",__LINE__,#c);return 1;}}while(0)
static std::vector<double> project(int na,int nb,int np,int*nc,int**col,double**J){
 int stride=na*nb,N=stride*np;std::vector<double>A(size_t(N)*N,0);
 for(int k=0;k<np;k++)for(int p=0;p<stride;p++){
  int row=p+stride*k;
  for(int q=0;q<nc[row];q++){int c=col[row][q]%stride,shift=col[row][q]/stride-k;double value=J[row][q]/(2*np);
   for(int t=0;t<np;t++){int plus=(t+shift+np)%np,minus=(t-shift+np)%np;A[size_t(p+stride*t)*N+c+stride*plus]+=value;A[size_t(p+stride*t)*N+c+stride*minus]+=value;}
  }
 }return A;
}
int main(){
 initialize();TwoPunctures_params_set_default();params_set_int("verbose",0);params_set_real("par_b",3);params_set_real("par_m_plus",.6);params_set_real("par_m_minus",.4);
 params_set_real("par_P_plus2",.036);params_set_real("par_S_plus3",.108);params_set_real("par_S_minus1",-.032);
 const int shape[3]={6,8,8},na=shape[0],nb=shape[1],np=shape[2],N=na*nb*np;
 derivs*v,*u,*dv;allocate_derivs(&v,N);allocate_derivs(&u,N);allocate_derivs(&dv,N);
 std::vector<double>F(N),native(N),answer(N),known(N),rhs(N),coeff(10*N,0);
 int*nc=ivector(0,N-1),**cols=imatrix(0,N-1,0,StencilSize-1);double**J=dmatrix(0,N-1,0,StencilSize-1);
 for(int p=0;p<N;p++)v->d0[p]=1e-4*std::cos(.19*p);
 // Full residual side effects, including raw AB/phi and Cartesian jets.
 derivs*kv,*ku;allocate_derivs(&kv,N);allocate_derivs(&ku,N);
 auto*workspace=Puncture_kokkos_BY_create(na,nb,np);CHECK(workspace);
 derivs malformed{};malformed.size=N;CHECK(Puncture_kokkos_BY_residual(workspace,&malformed,answer.data(),ku)==PK_INVALID);
 CHECK(Puncture_kokkos_BY_residual(workspace,kv,answer.data(),&malformed)==PK_INVALID);
 CHECK(Puncture_kokkos_BY_residual(workspace,kv,answer.data(),kv)==PK_INVALID);
 derivs alias=*ku;alias.d1=kv->d1;CHECK(Puncture_kokkos_BY_residual(workspace,kv,answer.data(),&alias)==PK_INVALID);
 alias=*ku;alias.d1=kv->d1+1;CHECK(Puncture_kokkos_BY_residual(workspace,kv,answer.data(),&alias)==PK_INVALID);
 CHECK(Puncture_kokkos_BY_residual(workspace,kv,ku->d0,ku)==PK_INVALID);
 for(int trial=0;trial<4;trial++){
  if(trial==3){params_set_real("par_m_plus",.55);params_set_real("par_m_minus",.45);}
  for(int p=0;p<N;p++)v->d0[p]=kv->d0[p]=trial==0?0:trial==1?1e-4*std::cos(.19*p):1e-4*(std::sin(1.234*p)+std::cos(Pi*(p/(na*nb))));
  F_of_v(1,na,nb,np,v,F.data(),u);CHECK(!Puncture_kokkos_BY_residual(workspace,kv,answer.data(),ku));
  double*ra[10]={v->d0,v->d1,v->d2,v->d3,v->d11,v->d12,v->d13,v->d22,v->d23,v->d33};
  double*rb[10]={kv->d0,kv->d1,kv->d2,kv->d3,kv->d11,kv->d12,kv->d13,kv->d22,kv->d23,kv->d33};
  double*ca[10]={u->d0,u->d1,u->d2,u->d3,u->d11,u->d12,u->d13,u->d22,u->d23,u->d33};
  double*cb[10]={ku->d0,ku->d1,ku->d2,ku->d3,ku->d11,ku->d12,ku->d13,ku->d22,ku->d23,ku->d33};
  for(int p=0;p<N;p++){CHECK(std::abs(F[p]-answer[p])/(1+std::abs(F[p]))<=1e-10);for(int d=0;d<10;d++){CHECK(std::abs(ra[d][p]-rb[d][p])/(1+std::abs(ra[d][p]))<=1e-10);CHECK(std::abs(ca[d][p]-cb[d][p])/(1+std::abs(ca[d][p]))<=1e-10);}}
 }
 params_set_real("par_b",4);CHECK(Puncture_kokkos_BY_residual(workspace,kv,answer.data(),ku)==PK_INVALID);params_set_real("par_b",3);
 params_set_real("par_P_plus2",.04);CHECK(Puncture_kokkos_BY_residual(workspace,kv,answer.data(),ku)==PK_INVALID);params_set_real("par_P_plus2",.036);
 params_set_real("par_S_plus3",.11);CHECK(Puncture_kokkos_BY_residual(workspace,kv,answer.data(),ku)==PK_INVALID);params_set_real("par_S_plus3",.108);
 params_set_real("par_m_plus",.6);params_set_real("par_m_minus",.4);
 // Independent original JVP versus a centered difference of the new F.
 for(int p=0;p<N;p++){v->d0[p]=1e-4*std::cos(.19*p);dv->d0[p]=std::sin(.31*p);}
 F_of_v(1,na,nb,np,v,F.data(),u);J_times_dv(1,na,nb,np,dv,native.data(),u);
 const double epsilon=1e-6;
 for(int p=0;p<N;p++)kv->d0[p]=v->d0[p]+epsilon*dv->d0[p];CHECK(!Puncture_kokkos_BY_residual(workspace,kv,answer.data(),ku));
 for(int p=0;p<N;p++)kv->d0[p]=v->d0[p]-epsilon*dv->d0[p];CHECK(!Puncture_kokkos_BY_residual(workspace,kv,rhs.data(),ku));
 for(int p=0;p<N;p++)CHECK(std::abs((answer[p]-rhs[p])/(2*epsilon)-native[p])/(1+std::abs(native[p]))<=1e-7);
 Puncture_kokkos_BY_destroy(workspace);free_derivs(kv);free_derivs(ku);
 for(int p=0;p<N;p++)v->d0[p]=1e-4*std::cos(.19*p);
 F_of_v(1,na,nb,np,v,F.data(),u);SetMatrix_JFD(1,na,nb,np,u,nc,cols,J);
 TP_Modal*m=TP_modal_create(1,na,nb,np,nc,cols,J);CHECK(m);auto M=import_BY_modal(m);auto projected=project(na,nb,np,nc,cols,J);
 const double*lu,*transfer,*lower,*scale,*fourier;const size_t*perm;CHECK(!TP_modal_factors(m,&lu,&transfer,&lower,&scale,&fourier,&perm));
 size_t nonidentity=0;for(int p=0;p<(np/2+1)*nb*na;p++)nonidentity+=perm[p]!=size_t(p%na);
 Operator A(shape,1,false);initialize_spectral(A,shape,3.);
 for(int k=0;k<np;k++)for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
  int p=i+na*(j+nb*k);double al=Pi*(i+.5)/na,be=Pi*(j+.5)/nb,a=.5*(1-std::cos(al)),B=-std::cos(be),X=2*std::atanh(a),R=Pih+2*std::atan(B),phi=2*Pi*k/np;
  double x=3*std::cosh(X)*std::cos(R),r=3*std::sinh(X)*std::sin(R),y=r*std::cos(phi),z=r*std::sin(phi),rp=std::sqrt((x-3)*(x-3)+y*y+z*z),rm=std::sqrt((x+3)*(x+3)+y*y+z*z);
  double psi=1+.3/rp+.2/rm+u->d0[p],p2=psi*psi,p4=p2*p2,w=std::pow(std::sin(al)*std::sin(be),3);
  coeff[10*p]=-.875*BY_KKofxyz(x,y,z)/(p4*p4)*w;coeff[10*p+4]=coeff[10*p+7]=coeff[10*p+9]=w;
 }A.coefficients=upload(coeff,"BY test coefficients");View output("BY test result",N);
 double largest_inverse=0,largest_action=0,largest_jvp=0;
 for(int mode=0;mode<np;mode++){
  for(int k=0;k<np;k++)for(int p=0;p<na*nb;p++)known[p+na*nb*k]=(1+.1*std::cos(.13*p))*std::cos(2*Pi*mode*k/np+.31);
  for(int row=0;row<N;row++){rhs[row]=0;for(int col=0;col<N;col++)rhs[row]+=projected[size_t(row)*N+col]*known[col];}
  auto in=upload(rhs,"BY manufactured RHS");M.apply(in,output);download(output,answer.data());double residual=0,rscale=0,error=0;
  for(int row=0;row<N;row++){double value=0;for(int col=0;col<N;col++)value+=projected[size_t(row)*N+col]*answer[col];residual=std::max(residual,std::abs(value-rhs[row]));rscale=std::max(rscale,std::abs(rhs[row]));error=std::max(error,std::abs(answer[row]-known[row]));}
  largest_inverse=std::max(largest_inverse,error);largest_action=std::max(largest_action,residual/std::max(1.,rscale));CHECK(error<=2e-9);CHECK(residual/std::max(1.,rscale)<=1e-10);
  std::copy(known.begin(),known.end(),dv->d0);J_times_dv(1,na,nb,np,dv,native.data(),u);in=upload(known,"BY JVP direction");A.apply(in,output);download(output,answer.data());
  for(int p=0;p<N;p++){double error=std::abs(native[p]-answer[p])/(1+std::abs(native[p]));largest_jvp=std::max(largest_jvp,error);CHECK(error<=1e-10);}
 }
 // The production BY bridge accepts both policies and their inherited host
 // callbacks. Manufacture an exact solution using the original spectral J.
 std::copy(known.begin(),known.end(),dv->d0);J_times_dv(1,na,nb,np,dv,rhs.data(),u);
 for(int method:{PK_GMRES,PK_BICGSTAB}){
  PK_Options options{method,300,32,1,0,0,1e-8,nullptr,nullptr};PK_Result result{};std::fill(answer.begin(),answer.end(),0.);
  if(method==PK_BICGSTAB){options.dot=forbidden_dot;options.norm=forbidden_norm;}
  CHECK(Puncture_kokkos_BY(na,nb,np,u->d0,m,rhs.data(),answer.data(),&options,&result,nullptr,nullptr)==PK_SUCCESS);
  std::copy(answer.begin(),answer.end(),dv->d0);J_times_dv(1,na,nb,np,dv,native.data(),u);double error=0;for(int p=0;p<N;p++)error+=(native[p]-rhs[p])*(native[p]-rhs[p]);CHECK(std::sqrt(error)<=1.1e-8);
 }
 // A separate well-conditioned cyclic block system forces nontrivial
 // pivots and dense lower blocks; the ordinary BY fixture needed no pivot.
 for(int k=0;k<np;k++)for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
  int row=i+na*(j+nb*k);nc[row]=0;
  auto append=[&](int col,double value){if(nc[row]>=StencilSize)std::abort();int q=nc[row]++;cols[row][q]=col;J[row][q]=value;};
  for(int q=0;q<na;q++){
   double value=i==q?4:.03*std::cos(i+q+.1*j);
   if(i<2&&q<2)value=i==q?.01:4;
   append(q+na*(j+nb*k),value);
   if(j)append(q+na*(j-1+nb*k),.02*std::sin(i+q+.2));
   if(j+1<nb&&q==i)append(q+na*(j+1+nb*k),.015*std::cos(i+q+.3));
  }
  append(i+na*(j+nb*((k+1)%np)),.01);append(i+na*(j+nb*((k+np-1)%np)),.01);CHECK(nc[row]<=StencilSize);
 }
 TP_Modal*pivoted=TP_modal_create(1,na,nb,np,nc,cols,J);CHECK(pivoted);auto pivotedM=import_BY_modal(pivoted);auto pivotedA=project(na,nb,np,nc,cols,J);
 CHECK(!TP_modal_factors(pivoted,&lu,&transfer,&lower,&scale,&fourier,&perm));nonidentity=0;for(int p=0;p<(np/2+1)*nb*na;p++)nonidentity+=perm[p]!=size_t(p%na);CHECK(nonidentity>0);
 for(int mode=0;mode<np;mode++){
  for(int k=0;k<np;k++)for(int p=0;p<na*nb;p++)known[p+na*nb*k]=(1+.1*std::cos(.13*p))*std::cos(2*Pi*mode*k/np+.31);
  for(int row=0;row<N;row++){rhs[row]=0;for(int col=0;col<N;col++)rhs[row]+=pivotedA[size_t(row)*N+col]*known[col];}
  auto in=upload(rhs,"pivoted full-lower RHS");pivotedM.apply(in,output);download(output,answer.data());
  for(int row=0;row<N;row++){CHECK(std::abs(answer[row]-known[row])<=1e-10);double value=0;for(int col=0;col<N;col++)value+=pivotedA[size_t(row)*N+col]*answer[col];CHECK(std::abs(value-rhs[row])/(1+std::abs(rhs[row]))<=1e-10);}
 }
 TP_modal_destroy(pivoted);
 TP_modal_destroy(m);free_derivs(v);free_derivs(u);free_derivs(dv);free_ivector(nc,0,N-1);free_imatrix(cols,0,N-1,0,StencilSize-1);free_dmatrix(J,0,N-1,0,StencilSize-1);
 printf("BY %s actual kernels: %d checks, M solution %.3e M action %.3e JVP %.3e\n",Exec::name(),checks,largest_inverse,largest_action,largest_jvp);return 0;
}
