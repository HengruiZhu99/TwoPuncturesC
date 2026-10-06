#ifndef PUNCTURE_KOKKOS_HPP
#define PUNCTURE_KOKKOS_HPP
#include <Kokkos_Core.hpp>
#include "PunctureExecution.h"
#include "TP_Modal.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>
#include <chrono>
#include <mutex>

namespace hispid {struct AxisDerivatives;}
namespace puncture {
using Exec=Kokkos::DefaultExecutionSpace;
using View=Kokkos::View<double*,Exec>;
using Indices=Kokkos::View<int*,Exec>;
using Host=Kokkos::View<const double*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using HostOut=Kokkos::View<double*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using Range=Kokkos::RangePolicy<Exec>;
using Action=std::function<void(View,View)>;
PunctureExecutionStats& statistics();
std::recursive_mutex& execution_mutex();
using ExecutionLock=std::unique_lock<std::recursive_mutex>;
void initialize();
inline double seconds(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
struct Timed {
 double&result,start;
 explicit Timed(double&v):result(v){Exec().fence();start=seconds();}
 ~Timed(){Exec().fence();result+=seconds()-start;}
};
struct TimedIf {
 double*result,start=0;
 TimedIf(double&v,bool enabled):result(enabled?&v:nullptr){if(result){Exec().fence();start=seconds();}}
 ~TimedIf(){if(result){Exec().fence();*result+=seconds()-start;}}
};
inline View upload(const double*v,size_t n,const char*label){
 View out(Kokkos::view_alloc(Kokkos::WithoutInitializing,std::string(label)),n);
 Timed t(statistics().transfer_seconds);Kokkos::deep_copy(out,Host(v,n));
 statistics().host_to_device_bytes+=n*sizeof(double);return out;
}
inline View upload(const std::vector<double>&v,const char*label){return upload(v.data(),v.size(),label);}
inline void download(View v,double*out){
 Timed t(statistics().transfer_seconds);Kokkos::deep_copy(HostOut(out,v.extent(0)),v);
 statistics().device_to_host_bytes+=v.extent(0)*sizeof(double);
}
inline void copy(View out,View in){Kokkos::deep_copy(out,in);}
inline void zero(View out){Kokkos::deep_copy(out,0.);}
template<class F> void each(View v,F f){Kokkos::parallel_for("vector",Range(0,v.extent(0)),f);}
inline double dot(View a,View b){double sum=0;Kokkos::parallel_reduce("dot",Range(0,a.extent(0)),KOKKOS_LAMBDA(int i,double&v){v+=a(i)*b(i);},sum);return sum;}
inline double length(View a){return std::sqrt(dot(a,a));}
inline bool finite(View a){int bad=0;Kokkos::parallel_reduce("finite",Range(0,a.extent(0)),KOKKOS_LAMBDA(int i,int&v){v+=!Kokkos::isfinite(a(i));},bad);return bad==0;}
inline void axpy(View out,double a,View in){each(out,KOKKOS_LAMBDA(int i){out(i)+=a*in(i);});}
inline void scale_copy(View out,double a,View in){each(out,KOKKOS_LAMBDA(int i){out(i)=a*in(i);});}
inline int action(const Action&f,View in,View out,int&count,double&timer){
 ++count;Timed t(timer);try{f(in,out);}catch(...){return PK_CALLBACK;}return finite(out)?PK_SUCCESS:PK_NONFINITE;
}
inline int actual(View b,View x,View ax,View r,const Action&A,PK_Result&result){
 if(!finite(x))return PK_NONFINITE;
 int status=action(A,x,ax,result.operator_calls,statistics().operator_seconds);if(status)return status;
 each(r,KOKKOS_LAMBDA(int i){r(i)=b(i)-ax(i);});
 result.true_residual=length(r);double rhs=length(b);result.relative_residual=rhs?result.true_residual/rhs:result.true_residual;
 return std::isfinite(result.true_residual)?PK_SUCCESS:PK_NONFINITE;
}
/* The C and Kokkos policies instantiate the same controller. */
int solve(View b,View x,const PK_Options&,const Action&A,const Action&M,PK_Result&,PK_Monitor=nullptr,void*monitor_context=nullptr);

/* Tensor-product kernels: contiguous output threads, no nested launches,
 * preallocated scratch, constant-annihilating off-diagonal differences. */
struct Spectral {
 int na,nb,np,nv,total;bool regular;
 View D[3],D2[3],inverse,inverse_phi,inverse_phi2,coordinate[2],coefficient[2],mapping;
 View work[10],scratch[6];
 Spectral(const int*n,int components,bool modal):na(n[0]),nb(n[1]),np(n[2]),nv(components),total(na*nb*np*nv),regular(modal){
  for(auto&w:work)w=View(Kokkos::view_alloc(Kokkos::WithoutInitializing,"derivative"),total);
  if(regular)for(auto&w:scratch)w=View(Kokkos::view_alloc(Kokkos::WithoutInitializing,"spectral scratch"),total);
 }
 void along(int axis,View matrix,View in,View out,bool difference){
  const int stride=axis==0?1:axis==1?na:na*nb,N=axis==0?na:axis==1?nb:np,V=nv;
  Kokkos::parallel_for("spectral line",Range(0,total),KOKKOS_LAMBDA(int q){
   int p=q/V,v=q%V,i=(p/stride)%N,base=p-i*stride;double value=0,center=in(q);
   for(int j=0;j<N;j++)if(!difference||j!=i)value+=matrix(i*N+j)*(in(V*(base+j*stride)+v)-(difference?center:0));out(q)=value;
  });
 }
 void apply(View input){
  copy(work[0],input);
  if(!regular){
   along(0,D[0],input,work[1],true);along(1,D[1],input,work[2],true);along(2,D[2],input,work[3],true);
   along(0,D2[0],input,work[4],true);along(1,D2[1],input,work[7],true);along(2,D2[2],input,work[9],true);
   along(1,D[1],work[1],work[5],true);along(2,D[2],work[1],work[6],true);along(2,D[2],work[2],work[8],true);return;
  }
  copy(scratch[0],input);along(0,D[0],input,scratch[1],true);along(1,D[1],input,scratch[2],true);
  along(0,D2[0],input,scratch[3],true);along(1,D2[1],input,scratch[4],true);along(0,D[0],scratch[2],scratch[5],true);
  auto P=scratch[0],Pa=scratch[1],Pb=scratch[2],Paa=scratch[3],Pbb=scratch[4],Pab=scratch[5],map=mapping;
  auto u=work[0],a=work[1],b=work[2],aa=work[4],bb=work[7],ab=work[5];const int V=nv;
  Kokkos::parallel_for("regular prolate modes",Range(0,total),KOKKOS_LAMBDA(int q){int p=q/V;const int z=15*p;
   u(q)=map(z)*P(q);a(q)=map(z+1)*Pa(q)+map(z+2)*P(q);b(q)=map(z+3)*Pb(q)+map(z+4)*P(q);
   aa(q)=map(z+5)*Paa(q)+map(z+6)*Pa(q)+map(z+7)*P(q);bb(q)=map(z+8)*Pbb(q)+map(z+9)*Pb(q)+map(z+10)*P(q);
   ab(q)=map(z+11)*Pab(q)+map(z+12)*Pb(q)+map(z+13)*Pa(q)+map(z+14)*P(q);
  });
  along(2,inverse_phi,work[0],work[3],false);along(2,inverse_phi2,work[0],work[9],false);
  along(2,inverse_phi,work[1],work[6],false);along(2,inverse_phi,work[2],work[8],false);
  for(int d:{0,1,2,4,7,5}){along(2,inverse,work[d],scratch[0],false);copy(work[d],scratch[0]);}
 }
};

/* Frozen per-point Jacobian in Cartesian jets, shared across equation systems.
 * Transformation coefficients are compact meridional data, never a dense
 * global matrix. The same coordinate chain as TP_CoordTransf is used. */
KOKKOS_INLINE_FUNCTION void prolate_coordinate_cache(double A,double B,double b,double*c,double*xrho){
 const double a=.5*(A+1),X=2*std::atanh(a),R=.5*std::acos(-1.0)+2*std::atan(B);
 const double sx=std::sinh(X),cx=std::cosh(X),cr=std::cos(R),sr=std::sin(R);
 // Match TP_CoordTransf's complex hyperbolic products before scaling by b.
 const double x=(cx*cr)*b,rho=(sx*sr)*b,u=(sx*cr)*b,v=(cx*sr)*b;
 const double den=u*u+v*v,re=u/den,im=-v/den;
 const double square_re=re*re-im*im,square_im=re*im+im*re;
 const double product_re=re*x-im*rho,product_im=re*rho+im*x;
 c[0]=re;c[1]=im;c[2]=-(square_re*product_re-square_im*product_im);
 c[3]=-(square_re*product_im+square_im*product_re);c[4]=re*re+im*im;c[5]=1/rho;
 xrho[0]=x;xrho[1]=rho;
}
KOKKOS_INLINE_FUNCTION void coordinate_chain(double A,double B,const double*c,double co,double si,double*f){
 double a=.5*(A+1),ax=1-a*a,axx=-a*ax,br=.5*(1+B*B),brr=B*br;
 double old[10];for(int d=0;d<10;d++)old[d]=f[d];
 for(int d=0;d<10;d++)f[d]*=A-1;
 f[1]+=old[0];f[4]+=2*old[1];f[5]+=old[2];f[6]+=old[3];
 f[4]=ax*ax*f[4]+axx*f[1];f[5]*=ax*br;f[6]*=ax;f[7]=br*br*f[7]+brr*f[2];f[8]*=br;f[1]*=ax;f[2]*=br;
 double cr=c[0],ci=c[1],ccr=c[2],cci=c[3],abs2=c[4],ri=c[5],ri2=ri*ri;
 double X=.5*f[1],Y=-.5*f[2],XX=.25*(f[4]-f[7]),YY=-.5*f[5],cb=.25*(f[4]+f[7])*abs2;
 double c2r=cr*cr-ci*ci,c2i=2*cr*ci,qr=ccr*X-cci*Y+c2r*XX-c2i*YY,qi=ccr*Y+cci*X+c2r*YY+c2i*XX;
 double ux=2*(X*cr-Y*ci),ur=-2*(X*ci+Y*cr),uxx=2*(cb+qr),urr=2*(cb-qr),uxr=-2*qi;
 double ux3=f[6]*cr+f[8]*ci,ur3=-f[6]*ci+f[8]*cr,u3=f[3],u33=f[9],cs2=co*co,sn2=si*si,s2=2*co*si,c2=cs2-sn2;
 f[1]=ux;f[2]=ur*co-u3*ri*si;f[3]=ur*si+u3*ri*co;f[4]=uxx;
 f[5]=uxr*co-ux3*ri*si;f[6]=uxr*si+ux3*ri*co;
 f[7]=urr*cs2+ri2*sn2*(u33+ur/ri)+s2*ri2*(u3-ur3/ri);
 f[8]=.5*s2*(urr-ri*ur-ri2*u33)-c2*ri2*(u3-ur3/ri);
 f[9]=urr*sn2+ri2*cs2*(u33+ur/ri)-s2*ri2*(u3-ur3/ri);
}
struct Operator {
 Spectral spectral;View coefficients,chain,trig,positions,fields;
 Operator(const int*n,int nv,bool regular):spectral(n,nv,regular),fields(Kokkos::view_alloc(Kokkos::WithoutInitializing,"Cartesian fields"),size_t(n[0])*n[1]*n[2]*nv*10){}
 void physical_fields(View in){
  spectral.apply(in);auto s=spectral;auto coef=coefficients,c=chain,t=trig;const int points=s.total/s.nv,na=s.na,nb=s.nb,V=s.nv;
  auto w0=s.work[0],w1=s.work[1],w2=s.work[2],w3=s.work[3],w4=s.work[4],w5=s.work[5],w6=s.work[6],w7=s.work[7],w8=s.work[8],w9=s.work[9];auto ca=s.coordinate[0],cb=s.coordinate[1],fout=fields;
  Kokkos::parallel_for("Cartesian derivatives",Range(0,points),KOKKOS_LAMBDA(int p){
   int i=p%na,j=(p/na)%nb,k=p/(na*nb);double f[4][10];
   for(int v=0;v<V;v++){int q=V*p+v;f[v][0]=w0(q);f[v][1]=w1(q);f[v][2]=w2(q);f[v][3]=w3(q);f[v][4]=w4(q);f[v][5]=w5(q);f[v][6]=w6(q);f[v][7]=w7(q);f[v][8]=w8(q);f[v][9]=w9(q);coordinate_chain(ca(i),cb(j),c.data()+6*(i+na*j),t(2*k),t(2*k+1),f[v]);}
   for(int v=0;v<V;v++)for(int d=0;d<10;d++)fout((p*V+v)*10+d)=f[v][d];
  });
 }
 void apply(View in,View out){
  physical_fields(in);auto f=fields,c=coefficients;const int V=spectral.nv,points=spectral.total/V;
  Kokkos::parallel_for("coupled Jacobian",Range(0,points),KOKKOS_LAMBDA(int p){for(int e=0;e<V;e++){double value=0;for(int v=0;v<V;v++)for(int d=0;d<10;d++)value+=c(((p*V+e)*V+v)*10+d)*f((p*V+v)*10+d);out(V*p+e)=value;}});
 }
};
void initialize_spectral(Operator&,const int*,double b,const std::vector<double>*coordinates=nullptr);
/* Builds mapped nodes, derivative/transform matrices and coordinate chains
 * in Exec. Host copies are only the small tables used by sampling/modal LU. */
void initialize_hispid_spectral(Operator&,hispid::AxisDerivatives&,const int*,double b);

/* Block inverse data is shared by Fourier partners/components. On CUDA,
 * invert the small pivoted LU blocks once, then use row-parallel GEMV in
 * the polar recurrence. OpenMP uses direct triangular solves per mode. */
struct ModalPointers {const double*lu,*transfer,*lower;const size_t*permutation;};
struct Modal {
 int na,nb,np,nv,groups;bool modal_output,diagonal_lower;
 bool compensated_projection=false; // opt-in by HiSpID; BY defaults unchanged
 View lu,transfer,lower,row_scale,forward,workspace,column,inverse;
 Indices permutation,block;
 Kokkos::View<ModalPointers*,Exec>native_blocks; // borrowed only during CPU apply
 Modal(int a,int b,int p,int v,int g,bool output,bool diag):na(a),nb(b),np(p),nv(v),groups(g),modal_output(output),diagonal_lower(diag){
  workspace=View("modal workspace",a*b*p*v);column=View("modal RHS",a*b*p*v);
 }
 void prepare();
 // projected=true accepts modal RHS values after Fourier projection and row
 // scaling, enabling inverse-only diagnostics without duplicating the kernel.
 void apply(View in,View out,bool projected=false);
};
Modal import_BY_modal(TP_Modal*);
}
#endif
