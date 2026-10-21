#include "HiSpID_internal.hpp"
#include "HiSpID_free_data.hpp"
#include "HiSpID_cache_kernels.hpp"
#include "HiSpID_axis.hpp"
#include "HiSpID_galerkin.hpp"
#include "PuncturePolynomialPreconditioner.hpp"
#include "PunctureLeftGMRES.hpp"
#ifndef HISPID_LEFT_GMRES
#define HISPID_LEFT_GMRES 0
#endif
#ifndef HISPID_DAMPED_MODAL_RELAXATION
#define HISPID_DAMPED_MODAL_RELAXATION 0
#endif
#include <memory>
#ifndef HISPID_OVERINTEGRATED_GALERKIN
#define HISPID_OVERINTEGRATED_GALERKIN 0
#endif
#include "HiSpID_tau.hpp"
#include "HiSpID_symmetry.hpp"
#include "HiSpID_modal_block.hpp"
#include "PunctureKrylov.h"
#include "PunctureExecution.h"
#include "PunctureModalProjection.hpp"
#ifndef HISPID_CONSTRAINT_BLOCK_PRECONDITIONER
#define HISPID_CONSTRAINT_BLOCK_PRECONDITIONER 0
#endif
#ifndef HISPID_COMPENSATED_MODAL_PROJECTION
#define HISPID_COMPENSATED_MODAL_PROJECTION 0
#endif
#ifdef PUNCTURES_KOKKOS
#include "PunctureKokkos.hpp"
#include <memory>
#include <mutex>
#define HI_INLINE KOKKOS_INLINE_FUNCTION
#else
#define HI_INLINE
#endif
#include <algorithm>
#include <chrono>
#include <cstring>
#include <map>
#include <numeric>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_linalg.h>
#ifndef HISPID_HALF_TURN_SOLVE
#define HISPID_HALF_TURN_SOLVE 0
#endif
#ifndef HISPID_ROW_POWER
#define HISPID_ROW_POWER 6
#endif
static_assert(HISPID_ROW_POWER==3||HISPID_ROW_POWER==6,"row power must be3or6");
#ifndef HISPID_MONOTONE_PRECONDITIONER
#define HISPID_MONOTONE_PRECONDITIONER 0
#endif
#ifndef HISPID_SPECTRAL_RADIAL_PRECONDITIONER
#define HISPID_SPECTRAL_RADIAL_PRECONDITIONER 0
#endif
#include "HiSpID_radial_preconditioner.hpp"
#include "HiSpID_curved_preconditioner.hpp"
#ifndef HISPID_CURVED_MODAL_PRECONDITIONER
#define HISPID_CURVED_MODAL_PRECONDITIONER 0
#endif
#ifndef HISPID_NEWTON_BACKTRACKS
#define HISPID_NEWTON_BACKTRACKS 10
#endif
static_assert(HISPID_NEWTON_BACKTRACKS>=10&&HISPID_NEWTON_BACKTRACKS<=40,
              "Newton backtracking depth must be between 10 and 40");
#ifndef HISPID_INFINITY_EQUILIBRATION
#define HISPID_INFINITY_EQUILIBRATION 0
#endif
static_assert(HISPID_INFINITY_EQUILIBRATION==0||HISPID_INFINITY_EQUILIBRATION==1,"infinity equilibration must be0or1");
extern "C" {
#include "TwoPunctures.h"
}
using hispid::Jet;
namespace {
const int hi[6]={1,1,1,2,2,3},hj[6]={1,2,3,2,3,3};
#ifdef PUNCTURES_KOKKOS
using Fields=Kokkos::Array<Kokkos::Array<double,10>,4>;
#else
using Fields=std::array<std::array<double,10>,4>;
#endif
using Cached=hispid::Cached;
using hispid::cache;
using hispid::L_and_div;
HI_INLINE double contraction(const Cached&c,const double*A,const double*B){
 double q=0;for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)for(int l=0;l<3;l++)
  q+=c.inv[3*i+k]*c.inv[3*j+l]*A[3*i+j]*B[3*k+l];
 return q;
}
HI_INLINE void eval(const Cached&c,const Fields&u,double*out,const Fields*du=nullptr){
 double psi=c.psi+u[0][0],A[9];if(!(psi>0)){for(int k=0;k<4;k++)out[k]=NAN;return;}
 for(int a=0;a<9;a++){
  A[a]=c.M[a];for(int k=0;k<3;k++)for(int d=0;d<4;d++)A[a]+=c.L[a][k][d]*u[k+1][d];
 }
 const Fields&v=du?*du:u;for(int k=0;k<4;k++)out[k]=0;
 for(int d=0;d<10;d++)out[0]+=c.lap[d]*v[0][d];
 for(int i=0;i<3;i++)for(int k=0;k<3;k++)for(int d=0;d<10;d++)out[i+1]+=c.vec[i][k][d]*v[k+1][d];
 double A2=contraction(c,A,A);
 if(!du){
#if HISPID_STABLE_SCALAR_SOURCE
  // Exact differences of powers avoid subtracting two background-sized terms.
  // L must be accumulated separately; A-M would lose small corrections.
  double L[9]={};for(int a=0;a<9;a++)for(int k=0;k<3;k++)for(int d=0;d<4;d++)L[a]+=c.L[a][k][d]*u[k+1][d];
  double sum5=0;for(int k=0;k<5;k++)sum5+=std::pow(psi,4-k)*std::pow(c.psi,k);
  const double ratio=c.psi/psi;double sum7=1,term=1;
  for(int k=1;k<7;k++){term*=ratio;sum7+=term;}
  const double inverse_difference=-u[0][0]/psi*sum7/std::pow(c.psi,7);
  out[0]+=c.g*(c.scalar_source-u[0][0]*c.R/8-c.K*c.K*u[0][0]*sum5/12
    +(2*contraction(c,c.M,L)+contraction(c,L,L))/(8*std::pow(psi,7))
    +c.seed_norm2*inverse_difference/8);
#else
  out[0]+=c.g*(-psi*c.R/8-std::pow(psi,5)*c.K*c.K/12+A2/(8*std::pow(psi,7))+c.lapPsi);
#endif
  for(int k=0;k<3;k++)out[k+1]+=c.g*(c.divM[k]-2.0/3*std::pow(psi,6)*c.gradK[k]);
 }else{
  double DA[9]={};for(int a=0;a<9;a++)for(int k=0;k<3;k++)for(int d=0;d<4;d++)DA[a]+=c.L[a][k][d]*v[k+1][d];
  out[0]+=c.g*((-c.R/8-5*std::pow(psi,4)*c.K*c.K/12-7*A2/(8*std::pow(psi,8)))*v[0][0]+contraction(c,A,DA)/(4*std::pow(psi,7)));
  for(int k=0;k<3;k++)out[k+1]-=c.g*4*std::pow(psi,5)*v[0][0]*c.gradK[k];
 }
}
#ifdef PUNCTURES_KOKKOS
struct HiKokkos {
 puncture::Operator op;
 Kokkos::View<Cached*,puncture::Exec>geometry;
 puncture::View base,input,output,averages;
#if HISPID_AXIS_TAU
 puncture::View tau_weight,tau_difference,focus_weight,focus_difference;
#endif
 bool raw_operator=false;
 double spectral_setup_seconds=0,geometry_setup_seconds=0;
 HiKokkos(const int*n,hispid::AxisDerivatives&axis,std::vector<Cached>&g,double b,const HiSpID_Config*execution_geometry=nullptr,int family=HISPID_SEED_QI,bool raw=false):
   op(n,4,true),
   base(Kokkos::view_alloc(Kokkos::WithoutInitializing,"frozen base fields"),40ull*n[0]*n[1]*n[2]),input("input",4ull*n[0]*n[1]*n[2]),output("output",4ull*n[0]*n[1]*n[2]),averages("azimuthal averages",2*n[0]*n[1]){
  raw_operator=raw;
#if HISPID_AXIS_TAU
  tau_weight=puncture::upload(hispid::tau_weights(n[0],n[1]),"axis endpoint weights");
  if(HISPID_FOCUS_TAU){focus_weight=puncture::upload(hispid::focus_bank(n[0],n[1],n[2]),"focus weights");focus_difference=puncture::View("focus differences",4*hispid::focus_rank*n[2]);}
  tau_difference=puncture::View("axis tau differences",4ull*hispid::tau_count(n[0],n[1])*n[2]);
#endif
  if(execution_geometry){initialize_execution_geometry(n,axis,b,*execution_geometry,family);return;}
  auto start=puncture::seconds();puncture::initialize_spectral(op,n,b,axis.coordinate.data());
  spectral_setup_seconds=puncture::seconds()-start;start=puncture::seconds();
  using H=Kokkos::View<const Cached*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  if constexpr(Kokkos::SpaceAccessibility<puncture::Exec,Kokkos::HostSpace>::accessible){
   // Borrow the context-owned cache on CPU: no second geometry allocation.
   geometry=decltype(geometry)(g.data(),g.size());
  }else{
   geometry=decltype(geometry)(Kokkos::view_alloc(Kokkos::WithoutInitializing,"background geometry"),g.size());
   puncture::Timed t(puncture::statistics().transfer_seconds);Kokkos::deep_copy(geometry,H(g.data(),g.size()));puncture::statistics().host_to_device_bytes+=g.size()*sizeof(Cached);
  }
  geometry_setup_seconds=puncture::seconds()-start;start=puncture::seconds();
  op.spectral.inverse=puncture::upload(axis.inverse,"inverse Fourier");op.spectral.inverse_phi=puncture::upload(axis.inverse_phi,"Fourier first");op.spectral.inverse_phi2=puncture::upload(axis.inverse_phi2,"Fourier second");
  std::vector<double>mapping(15*g.size());
  for(size_t p=0;p<g.size();p++){
   int i=p%n[0],j=(p/n[0])%n[1],mode=p/(n[0]*n[1]),m=mode<=n[2]/2?mode:mode-n[2]/2,r=hispid::AxisDerivatives::exponent(m);
   double a=.5*(axis.coordinate[0][i]+1),B=axis.coordinate[1][j],eta=-2*B/(1+B*B),den=1-eta*eta;
   double ma=axis.first_map[0][i],mma=axis.second_map[0][i],mb=axis.first_map[1][j],mmb=axis.second_map[1][j],S=hispid::AxisDerivatives::radial(a,r)*hispid::AxisDerivatives::angular(eta,r);
   double la=.5*(1/(1+a)+r/a),lla=-.25*(1/((1+a)*(1+a))+r/(a*a));
   double etaB=-2*(1-B*B)/std::pow(1+B*B,2),etaBB=4*B*(3-B*B)/std::pow(1+B*B,3),le=-r*eta/den,lle=-r*(1+eta*eta)/(den*den),lb=le*etaB,llb=lle*etaB*etaB+le*etaBB;
   double*q=mapping.data()+15*p;
   q[0]=S;q[1]=S*ma;q[2]=S*la;q[3]=S*mb;q[4]=S*lb;
   q[5]=S*ma*ma;q[6]=S*(mma+2*la*ma);q[7]=S*(la*la+lla);
   q[8]=S*mb*mb;q[9]=S*(mmb+2*lb*mb);q[10]=S*(lb*lb+llb);
   q[11]=S*ma*mb;q[12]=S*la*mb;q[13]=S*lb*ma;q[14]=S*la*lb;
  }
  op.spectral.mapping=puncture::upload(mapping,"regular factors");
  spectral_setup_seconds+=puncture::seconds()-start;
 }
 void initialize_execution_geometry(const int*n,hispid::AxisDerivatives&axis,double b,const HiSpID_Config&config,int family){
   auto start=puncture::seconds();puncture::initialize_hispid_spectral(op,axis,n,b);spectral_setup_seconds=puncture::seconds()-start;
   const int points=n[0]*n[1]*n[2],na=n[0],nb=n[1];
   geometry=decltype(geometry)(Kokkos::view_alloc(Kokkos::WithoutInitializing,"kernel background geometry"),points);
   auto geom=geometry;auto position=op.positions,trig=op.trig,ca=op.spectral.coordinate[0];puncture::Indices errors("geometry statuses",points);
   start=puncture::seconds();
   Kokkos::parallel_for("spinning boosted seed and operator cache",puncture::Range(0,points),KOKKOS_LAMBDA(int p){
    int row=p%(na*nb),i=row%na,j=row/na,k=p/(na*nb);double xyz[3]={position(2*row),position(2*row+1)*trig(2*k),position(2*row+1)*trig(2*k+1)};
    hispid::BackgroundT<double>bg;int status=hispid::background_geometry(config,xyz,bg,true,family);
    if(status){errors(p)=status;return;}
    cache(bg,geom(p));double sn=std::sin(Pih*(2*i+1)/na)*std::sin(Pih*(2*j+1)/nb);geom(p).weight=std::pow(sn,HISPID_ROW_POWER);
    const double a=.5*(ca(i)+1);if constexpr(HISPID_INFINITY_EQUILIBRATION){geom(p).weight/=std::pow(1-a*a,6);}
    errors(p)=hispid::finite_cache(geom(p))?hispid::geometry_ok:hispid::geometry_nonfinite;
   });
   puncture::Exec().fence();geometry_setup_seconds=puncture::seconds()-start;
   int failed=points;Kokkos::parallel_reduce("first failed geometry point",puncture::Range(0,points),KOKKOS_LAMBDA(int p,int&first){if(errors(p)&&p<first)first=p;},Kokkos::Min<int>(failed));
   if(failed<points){int status=0;Kokkos::deep_copy(status,Kokkos::subview(errors,failed));throw std::runtime_error("execution geometry failed at point "+std::to_string(failed)+" (status "+std::to_string(status)+")");}
 }
 void weight_rows(puncture::View out){
  auto c=geometry;Kokkos::parallel_for("Galerkin row equilibration",puncture::Range(0,geometry.extent(0)),KOKKOS_LAMBDA(int p){for(int v=0;v<4;v++)out(4*p+v)*=c(p).weight;});
 }
 void row_weights(std::vector<double>&out){
  const int rows=op.spectral.na*op.spectral.nb;puncture::View weights("row weight export",rows);auto c=geometry;
  Kokkos::parallel_for("row weight export",puncture::Range(0,rows),KOKKOS_LAMBDA(int p){weights(p)=c(p).weight;});out.resize(rows);puncture::download(weights,out.data());
 }
 void make_coefficients(const double*values,double*out){
  {puncture::Timed t(puncture::statistics().transfer_seconds);Kokkos::deep_copy(input,puncture::Host(values,input.extent(0)));puncture::statistics().host_to_device_bytes+=input.extent(0)*8;}
  op.spectral.along(0,op.spectral.coefficient[0],input,op.spectral.work[0],false);
  op.spectral.along(1,op.spectral.coefficient[1],op.spectral.work[0],output,false);puncture::download(output,out);
 }
 void fields(const double*v,bool reference){
  {puncture::Timed t(puncture::statistics().transfer_seconds);Kokkos::deep_copy(input,puncture::Host(v,input.extent(0)));puncture::statistics().host_to_device_bytes+=input.extent(0)*8;}
  fields(input,reference);
 }
 void fields(puncture::View values,bool reference){
  op.physical_fields(values);auto f=op.fields;auto c=geometry;auto out=base;
  Kokkos::parallel_for("reference correction",puncture::Range(0,geometry.extent(0)),KOKKOS_LAMBDA(int p){for(int v=0;v<4;v++)for(int d=0;d<10;d++)out(40*p+10*v+d)=f(40*p+10*v+d)+(reference&&v==0?c(p).far_correction[d]:0);});
 }
 void replace_axis_rows(puncture::View in,puncture::View out){
  if(raw_operator)return;
#if HISPID_AXIS_TAU
  const int na=op.spectral.na,nb=op.spectral.nb,np=op.spectral.np,count=hispid::tau_count(na,nb);
  auto w=tau_weight,d=tau_difference,inv=op.spectral.inverse;
  if(HISPID_FOCUS_TAU){
   auto fw=focus_weight,fd=focus_difference;const int N=na*nb;
   using Team=Kokkos::TeamPolicy<puncture::Exec>;using Member=Team::member_type;
   Kokkos::parallel_for("focus regularity moments",Team(4*hispid::focus_rank*np,puncture::reduction_team_size),KOKKOS_LAMBDA(const Member&team){
    int task=team.league_rank(),v=task%4,k=(task/4)%np,side=task/(4*np),m=k<=np/2?k:k-np/2;
    if(m>4){Kokkos::single(Kokkos::PerTeam(team),[&](){fd(task)=0;});return;}
    puncture::ProjectionAccumulator sum;Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,N),[&](int p,puncture::ProjectionAccumulator&q){q.product(fw((hispid::focus_rank*m+side)*N+p),in(4*(p+N*k)+v));},Kokkos::Sum<puncture::ProjectionAccumulator>(sum));
    Kokkos::single(Kokkos::PerTeam(team),[&](){int row=hispid::focus_row(side,na,nb);double old=0;for(int l=0;l<np;l++)old+=inv(l*np+k)*out(4*(row+N*l)+v);fd(task)=(sum.high+sum.low)-old;});
   });
  }
  Kokkos::parallel_for("axis tau differences",puncture::Range(0,4*count*np),KOKKOS_LAMBDA(int z){
   int v=z%4,m=(z/4)%np,q=z/(4*np);
   d(z)=hispid::tau_delta(q,m,v,na,nb,np,w.data(),inv.data(),in.data(),out.data());
  });
  Kokkos::parallel_for("axis tau rows",puncture::Range(0,4*count*np),KOKKOS_LAMBDA(int z){
   int v=z%4,k=(z/4)%np,q=z/(4*np),i,j;hispid::tau_node(q,na,nb,i,j);
   double correction=0;for(int m=0;m<np;m++)correction+=inv(k*np+m)*d(4*(q*np+m)+v);
   out(4*(i+na*(j+nb*k))+v)+=correction;
  });
  if(HISPID_FOCUS_TAU){auto fd=focus_difference;const int N=na*nb;
   Kokkos::parallel_for("focus regularity rows",puncture::Range(0,4*hispid::focus_rank*np),KOKKOS_LAMBDA(int task){
    int v=task%4,k=(task/4)%np,side=task/(4*np);double q=0;
    for(int m=0;m<np;m++)q+=inv(k*np+m)*fd(4*(side*np+m)+v);
    out(4*(hispid::focus_row(side,na,nb)+N*k)+v)+=q;
   });
  }
#endif
 }
 void residual(const double*v,double*r){
  {puncture::Timed t(puncture::statistics().transfer_seconds);Kokkos::deep_copy(input,puncture::Host(v,input.extent(0)));puncture::statistics().host_to_device_bytes+=input.extent(0)*8;}
  residual(input,output);puncture::download(output,r);
 }
 void residual(puncture::View in,puncture::View out){
  fields(in,true);auto f=base;auto c=geometry;const bool raw=raw_operator;
  Kokkos::parallel_for("HiSpID nonlinear equations",puncture::Range(0,geometry.extent(0)),KOKKOS_LAMBDA(int p){Fields u{};for(int v=0;v<4;v++)for(int d=0;d<10;d++)u[v][d]=f(40*p+10*v+d);double result[4];eval(c(p),u,result);for(int v=0;v<4;v++)out(4*p+v)=result[v]*(raw?1:c(p).weight);});
  replace_axis_rows(in,out);
 }
 void apply(puncture::View in,puncture::View out){
  op.physical_fields(in);auto f=op.fields,b=base;auto c=geometry;const bool raw=raw_operator;
  Kokkos::parallel_for("HiSpID Jacobian",puncture::Range(0,geometry.extent(0)),KOKKOS_LAMBDA(int p){Fields u{},du{};for(int v=0;v<4;v++)for(int d=0;d<10;d++){u[v][d]=b(40*p+10*v+d);du[v][d]=f(40*p+10*v+d);}double result[4];eval(c(p),u,result,&du);for(int v=0;v<4;v++)out(4*p+v)=result[v]*(raw?1:c(p).weight);});
  replace_axis_rows(in,out);
 }
 void azimuthal_averages(double*out){
  const int stride=op.spectral.na*op.spectral.nb,np=op.spectral.np;auto c=geometry;auto f=base,avg=averages;
  Kokkos::parallel_for("modal coefficient averages",puncture::Range(0,stride),KOKKOS_LAMBDA(int row){double mu=0,potential=0;
   for(int k=0;k<np;k++){int p=row+stride*k;const auto&g=c(p);double A[9];
    for(int a=0;a<9;a++){A[a]=g.M[a];for(int v=0;v<3;v++)for(int d=0;d<4;d++)A[a]+=g.L[a][v][d]*f(40*p+10*(v+1)+d);}
    double psi=g.psi+f(40*p);mu+=(g.inv[0]+g.inv[4]+g.inv[8])/(3*np);
    potential+=g.g*(-g.R/8-5*std::pow(psi,4)*g.K*g.K/12-7*contraction(g,A,A)/(8*std::pow(psi,8)))/np;
   }avg(2*row)=mu;avg(2*row+1)=potential;
  });puncture::download(averages,out);
 }
 void curved_symbols(double b,std::vector<hispid::ModalSymbol>&out){
  const int na=op.spectral.na,nb=op.spectral.nb,np=op.spectral.np,N=na*nb,groups=2*(np/2+1);
  Kokkos::View<hispid::ModalSymbol*,puncture::Exec> symbols("curved modal symbol",N*groups);
  auto c=geometry;auto ca=op.spectral.coordinate[0],cb=op.spectral.coordinate[1],trig=op.trig;
  Kokkos::parallel_for("curved modal averages",puncture::Range(0,N*groups),KOKKOS_LAMBDA(int q){
   int row=q%N,group=q/N;hispid::ModalSymbol avg;
   double a=.5*(ca(row%na)+1),t=a*a,C=(1+t)/(1-t),B=cb(row/na),eta=-2*B/(1+B*B),rho=b*std::sqrt((C*C-1)*(1-eta*eta));
   for(int k=0;k<np;k++){
    auto v=hispid::curved_symbol(c(row+N*k),b,b*C*eta,rho*trig(2*k),rho*trig(2*k+1),group/2,group%2);
    avg.zero+=v.zero/np;avg.mixed+=v.mixed/np;for(int a=0;a<2;a++){avg.first[a]+=v.first[a]/np;avg.second[a]+=v.second[a]/np;}
   }symbols(q)=avg;
  });
  auto host=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),symbols);out.assign(host.data(),host.data()+host.extent(0));
 }
 unsigned long long bytes()const{
  unsigned long long extra=0;
#if HISPID_AXIS_TAU
  extra=8*(tau_weight.extent(0)+tau_difference.extent(0)+focus_weight.extent(0)+focus_difference.extent(0));
#endif
  auto n=geometry.extent(0);return extra+n*(sizeof(Cached)+8*(40+40+8+64+15))+8*(op.chain.extent(0)+op.trig.extent(0)+op.positions.extent(0)+op.spectral.coefficient[0].extent(0)+op.spectral.coefficient[1].extent(0));
 }
};
// Rectangular tensor transfers remain on the selected execution backend.
struct DeviceGalerkinTransfer {
 puncture::View B[2],R[2],scratch;
 std::array<int,3>n,q;
 explicit DeviceGalerkinTransfer(const hispid::GalerkinTransfer&t):n(t.n),q(t.q){
  for(int a=0;a<2;a++){B[a]=puncture::upload(t.B[a],"Galerkin prolongation");R[a]=puncture::upload(t.R[a],"Galerkin projection");}
  scratch=puncture::View(Kokkos::view_alloc(Kokkos::WithoutInitializing,"Galerkin transfer scratch"),4ull*q[0]*q[1]*q[2]);
 }
 static void along(int axis,int rows,int cols,int radial,int polar,int phi,puncture::View matrix,puncture::View in,puncture::View out){
  const int ia=axis==0?cols:radial,ib=axis==0?polar:cols,oa=axis==0?rows:radial,ob=axis==0?polar:rows;
  Kokkos::parallel_for("Galerkin tensor transfer",puncture::Range(0,4*oa*ob*phi),KOKKOS_LAMBDA(int z){
   int v=z%4,p=z/4,i=p%oa,j=(p/oa)%ob,k=p/(oa*ob),row=axis==0?i:j;
   int anchor=Kokkos::min(cols-1,int((row+.5)*cols/rows));
   int base=axis==0?ia*(j+ib*k):i+ia*ib*k,stride=axis==0?1:ia;
   const double*x=in.data()+4*base+v;double center=x[4*stride*anchor];
   out(z)=center+puncture::compensated_projection(x,4*stride,matrix.data()+row*cols,cols,center);
  });
 }
 void prolong(puncture::View in,puncture::View out){
  along(0,q[0],n[0],n[0],n[1],n[2],B[0],in,scratch);
  along(1,q[1],n[1],q[0],n[1],n[2],B[1],scratch,out);
 }
 void restrict(puncture::View in,puncture::View out){
  along(1,n[1],q[1],q[0],q[1],q[2],R[1],in,scratch);
  along(0,n[0],q[0],q[0],n[1],q[2],R[0],scratch,out);
 }
 unsigned long long bytes()const{return 8ull*(scratch.extent(0)+B[0].extent(0)+B[1].extent(0)+R[0].extent(0)+R[1].extent(0));}
};

#endif
}
struct HiSpID_Delete { void operator()(HiSpID_Data*p)const{HiSpID_destroy(p);} };
struct HiSpID_Data {
 HiSpID_Config config,local;
 long double allocation_bound=0;
 int seed_family=HISPID_SEED_QI;
 HiSpID_FreeData free_data{};
 HiSpID_Wave22 wave{};
 bool has_wave=false;
 double origin[3],frame[3][3],b;
 int npt,ntotal;
 std::vector<Cached> geometry;
 std::vector<double> values,coefficients;
#if HISPID_AXIS_TAU
 std::vector<double> tau_weight,tau_difference,focus_weight,focus_difference;
#endif
 derivs *work=nullptr;
 hispid::AxisDerivatives derivatives;
 std::vector<Fields> basefields;
 std::vector<hispid::ModalSymbol>curved_symbols;
 HiSpID_Diagnostics diag{};
 bool coefficients_valid=false;
 bool sampler_only=false,operator_only=false;
 std::unique_ptr<HiSpID_Data,HiSpID_Delete>quadrature;
 std::unique_ptr<hispid::GalerkinTransfer>transfer;
 std::vector<double>quadrature_residual;
 bool axisymmetric=false,half_turn=false;
 HiSpID_SetupStatistics setup{int(sizeof(HiSpID_SetupStatistics)),0,std::numeric_limits<long double>::digits,0,0,0};
 int jvp_applications=0,preconditioner_applications=0;
 std::vector<std::array<double,4>>linear_history;
 double last_gmres_relative=0;
 HiSpID_SolveOptions resolved_options{};
#ifdef PUNCTURES_KOKKOS
 std::unique_ptr<HiKokkos>device;
 std::unique_ptr<DeviceGalerkinTransfer>device_transfer;
 std::vector<double>row_weights;
#endif
};
namespace {
#ifdef PUNCTURES_KOKKOS
puncture::ExecutionLock execution_lock(const HiSpID_Data*s){
 puncture::ExecutionLock lock(puncture::execution_mutex(),std::defer_lock);if(s&&s->device)lock.lock();return lock;
}
#endif
// Both cache construction and off-grid sampling use exactly the same free data.
void context_background(const HiSpID_Data&s,const double*x,hispid::Background&bg){
 hispid::background(s.local,x,bg,s.seed_family);
 if(!s.has_wave&&!s.free_data.evaluate)return;
 double lab[3];for(int i=0;i<3;i++){lab[i]=s.origin[i];for(int j=0;j<3;j++)lab[i]+=s.frame[i][j]*x[j];}
 hispid::Jet h[3][3],ht[3][3],local_h[3][3],local_ht[3][3];
 if(s.has_wave){
  int active=s.config.hole[0].mass>0?0:1;
  hispid::wave22(s.wave,s.config.hole[active].center,lab,h,ht);
 }else{
  double jets[2][90]={};
  if(s.free_data.evaluate(lab,jets[0],jets[1],s.free_data.user))throw std::runtime_error("free-data callback failed");
  constexpr int a[6]={0,0,0,1,1,2},b[6]={0,1,2,1,2,2};
  for(int t=0;t<2;t++)for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   const double*v=jets[t]+10*(3*i+j);hispid::Jet&q=t?ht[i][j]:h[i][j];
   for(int d=0;d<10;d++)if(!std::isfinite(v[d])||v[d]!=jets[t][10*(3*j+i)+d])throw std::runtime_error("nonfinite or asymmetric free-data tensor jets");
   q.v=v[0];for(int d=0;d<3;d++)q.d[d+1]=v[d+1];
   for(int d=0;d<6;d++)q.h[a[d]+1][b[d]+1]=q.h[b[d]+1][a[d]+1]=v[d+4];
  }
 }
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int a=0;a<3;a++)for(int b=0;b<3;b++){
  const long double tensor=s.frame[a][i]*s.frame[b][j];
  for(int t=0;t<2;t++){
   const hispid::Jet&q=t?ht[a][b]:h[a][b];hispid::Jet&v=t?local_ht[i][j]:local_h[i][j];
   v.v+=tensor*q.v;
   for(int d=0;d<3;d++)for(int c=0;c<3;c++){
    v.d[d+1]+=tensor*s.frame[c][d]*q.d[c+1];
    for(int e=0;e<3;e++)for(int f=0;f<3;f++)v.h[d+1][e+1]+=tensor*s.frame[c][d]*s.frame[f][e]*q.h[c+1][f+1];
   }
  }
 }
 hispid::perturb_background(s.local,s.seed_family,x,local_h,local_ht,bg);
}
int pindex(const HiSpID_Data&s,int i,int j,int k){return i+s.local.n[0]*(j+s.local.n[1]*k);}
void project_solve_sector(const HiSpID_Data&s,double*v,bool unknowns){
 // Only unknowns are restricted. Full RHS, JVP and nonlinear diagnostics
 // retain the original equations and original true-residual stopping norm.
 if(s.half_turn&&unknowns){
  const int N=s.local.n[0]*s.local.n[1],np=s.local.n[2];
  for(int q=0;q<s.ntotal;q++){const int slot=q/(4*N),mode=slot<=np/2?slot:slot-np/2;
   if(mode%2!=(q%4>=2?1:0))v[q]=0.;}
 }
 if(!s.axisymmetric)return;
 const int stride=s.local.n[0]*s.local.n[1],np=s.local.n[2];
 for(int row=0;row<stride;row++){
  if(unknowns)hispid::axisymmetric_unknown_row(v,row,stride,np);
  else hispid::axisymmetric_equation_row(v,row,stride,np);
 }
}
#ifdef PUNCTURES_KOKKOS
void project_solve_sector(const HiSpID_Data&s,puncture::View v,bool unknowns){
 if(s.half_turn&&unknowns){
  const int N=s.local.n[0]*s.local.n[1],np=s.local.n[2];
  puncture::each(v,KOKKOS_LAMBDA(int q){const int slot=q/(4*N),mode=slot<=np/2?slot:slot-np/2;
   if(mode%2!=(q%4>=2?1:0))v(q)=0.;});
 }
 if(!s.axisymmetric)return;
 const int stride=s.local.n[0]*s.local.n[1],np=s.local.n[2];
 Kokkos::parallel_for("axisymmetric sector",puncture::Range(0,stride),KOKKOS_LAMBDA(int row){
  if(unknowns)hispid::axisymmetric_unknown_row(v.data(),row,stride,np);
  else hispid::axisymmetric_equation_row(v.data(),row,stride,np);
 });
}
#endif
double row_weight(const HiSpID_Data&s,int p){
#ifdef PUNCTURES_KOKKOS
 if(s.device)return s.row_weights[p%(s.local.n[0]*s.local.n[1])];
#endif
 return s.geometry[p].weight;
}
double dot(const std::vector<double>&a,const std::vector<double>&b){double q=0;for(size_t i=0;i<a.size();i++)q+=a[i]*b[i];return q;}
double norm2v(const std::vector<double>&v){return std::sqrt(dot(v,v));}
double norminf(const std::vector<double>&v){double q=0;for(double a:v){if(!std::isfinite(a))return INFINITY;q=std::max(q,std::abs(a));}return q;}
void gather(derivs*w,int p,Fields&f){
 double *a[10]={w->d0,w->d1,w->d2,w->d3,w->d11,w->d12,w->d13,w->d22,w->d23,w->d33};
 for(int k=0;k<4;k++)for(int d=0;d<10;d++)f[k][d]=a[d][4*p+k];
}
void transform(HiSpID_Data&s,int i,int j,int k,Fields&f,double*x=nullptr){
 double A=s.derivatives.coordinate[0][i],B=s.derivatives.coordinate[1][j];
 double dat[10][4];derivs w{};double **ptr[10]={&w.d0,&w.d1,&w.d2,&w.d3,&w.d11,&w.d12,&w.d13,&w.d22,&w.d23,&w.d33};
 for(int d=0;d<10;d++){*ptr[d]=dat[d];for(int v=0;v<4;v++)dat[d][v]=(A-1)*f[v][d];}
 for(int v=0;v<4;v++){dat[1][v]+=f[v][0];dat[4][v]+=2*f[v][1];dat[5][v]+=f[v][2];dat[6][v]+=f[v][3];}
 double X,R,xx,rr,y,z;AB_To_XR(4,A,B,&X,&R,&w);C_To_c(4,X,R,&xx,&rr,s.b,&w);
 rx3_To_xyz(4,xx,rr,2*Pi*k/s.local.n[2],&y,&z,&w);
 for(int v=0;v<4;v++)for(int d=0;d<10;d++)f[v][d]=dat[d][v];
 if(x){x[0]=xx;x[1]=y;x[2]=z;}
}
void fields(HiSpID_Data&s,const double*values,std::vector<Fields>&f,bool include_reference=true){
#ifdef PUNCTURES_KOKKOS
 if(s.device){s.device->fields(values,include_reference);f.resize(s.npt);static_assert(sizeof(Fields)==40*sizeof(double));puncture::download(s.device->base,reinterpret_cast<double*>(f.data()));return;}
#endif
 std::copy(values,values+s.ntotal,s.work->d0);
 s.derivatives.apply(4,s.work);
 f.resize(s.npt);for(int k=0;k<s.local.n[2];k++)for(int j=0;j<s.local.n[1];j++)for(int i=0;i<s.local.n[0];i++){
  int p=pindex(s,i,j,k);gather(s.work,p,f[p]);transform(s,i,j,k,f[p]);
  if(include_reference)for(int d=0;d<10;d++)f[p][0][d]+=s.geometry[p].far_correction[d];
 }
}
void replace_axis_rows(HiSpID_Data&s,const double*v,double*r){
 if(s.operator_only)return;
#if HISPID_AXIS_TAU
 const auto*n=s.local.n;
 if(s.tau_weight.empty())s.tau_weight=hispid::tau_weights(n[0],n[1]);
 if(HISPID_FOCUS_TAU){
  if(s.focus_weight.empty())s.focus_weight=hispid::focus_bank(n[0],n[1],n[2]);
  hispid::focus_differences(n[0],n[1],n[2],s.focus_weight.data(),s.derivatives.inverse.data(),v,r,s.focus_difference);
 }
 hispid::tau_apply(n[0],n[1],n[2],s.tau_weight.data(),s.derivatives.inverse.data(),v,r,s.tau_difference);
 if(HISPID_FOCUS_TAU)hispid::focus_replace(n[0],n[1],n[2],s.derivatives.inverse.data(),s.focus_difference,r);
#endif
}
#ifdef PUNCTURES_KOKKOS
void device_apply(HiSpID_Data&s,puncture::View in,puncture::View out){
 if(s.quadrature){
  auto&q=*s.quadrature->device;s.device_transfer->prolong(in,q.input);
  q.apply(q.input,q.output);s.device_transfer->restrict(q.output,out);s.device->weight_rows(out);s.device->replace_axis_rows(in,out);
 }else s.device->apply(in,out);
}
#endif
void residual(HiSpID_Data&s,const double*v,double*r){
#ifdef PUNCTURES_KOKKOS
 if(s.device&&s.quadrature){
  auto&d=*s.device;auto&q=*s.quadrature->device;
  {puncture::Timed t(puncture::statistics().transfer_seconds);Kokkos::deep_copy(d.input,puncture::Host(v,s.ntotal));puncture::statistics().host_to_device_bytes+=8ull*s.ntotal;}
  s.device_transfer->prolong(d.input,q.input);q.residual(q.input,q.output);
  s.device_transfer->restrict(q.output,d.output);d.weight_rows(d.output);d.replace_axis_rows(d.input,d.output);puncture::download(d.output,r);return;
 }
#endif
 if(s.quadrature){
  auto&q=*s.quadrature;s.transfer->prolong(v,q.values.data());
  s.quadrature_residual.resize(q.ntotal);residual(q,q.values.data(),s.quadrature_residual.data());
  s.transfer->restrict(s.quadrature_residual.data(),r);for(int p=0;p<s.npt;p++)for(int k=0;k<4;k++)r[4*p+k]*=row_weight(s,p);replace_axis_rows(s,v,r);return;
 }
#ifdef PUNCTURES_KOKKOS
 if(s.device){s.device->residual(v,r);return;}
#endif
 fields(s,v,s.basefields);for(int p=0;p<s.npt;p++){
  eval(s.geometry[p],s.basefields[p],r+4*p);for(int k=0;k<4;k++)r[4*p+k]*=s.operator_only?1:s.geometry[p].weight;
 }
 replace_axis_rows(s,v,r);
}
void jvp(HiSpID_Data&s,const double*v,double*r){
 s.jvp_applications++;
#ifdef PUNCTURES_KOKKOS
 if(s.device){Kokkos::deep_copy(s.device->input,puncture::Host(v,s.ntotal));device_apply(s,s.device->input,s.device->output);puncture::download(s.device->output,r);return;}
#endif
 if(s.quadrature){
  auto&q=*s.quadrature;s.transfer->prolong(v,q.values.data());
  s.quadrature_residual.resize(q.ntotal);jvp(q,q.values.data(),s.quadrature_residual.data());
  s.transfer->restrict(s.quadrature_residual.data(),r);for(int p=0;p<s.npt;p++)for(int k=0;k<4;k++)r[4*p+k]*=row_weight(s,p);replace_axis_rows(s,v,r);return;
 }
#ifdef PUNCTURES_KOKKOS
 if(s.device){Kokkos::deep_copy(s.device->input,puncture::Host(v,s.ntotal));s.device->apply(s.device->input,s.device->output);puncture::download(s.device->output,r);return;}
#endif
 std::vector<Fields>d;fields(s,v,d,false);for(int p=0;p<s.npt;p++){
  eval(s.geometry[p],s.basefields[p],r+4*p,&d[p]);for(int k=0;k<4;k++)r[4*p+k]*=s.operator_only?1:s.geometry[p].weight;
 }
 replace_axis_rows(s,v,r);
}
/* Exact block elimination of the five-point modal FD approximation. Each
 * polar row is a dense radial block after elimination. Cosine/sine partners
 * share factors, as do the three approximate vector equations. This removes
 * the long-wavelength error of ILU(0); residuals/JVPs remain pseudospectral. */
using hispid::ModalBlock;
struct Sparse {
 std::vector<std::vector<int>>col;
 std::vector<std::vector<double>>val;
 std::vector<int>diag;
 const hispid::AxisDerivatives*modal=nullptr;
 std::vector<double>row_scale;
 std::vector<ModalBlock>blocks;
 int scalar_factorizations=0,vector_factorizations=0;
 // Axial unknowns use only m0 and m1. Higher output modes are removed by
 // project_solve_sector; alias their temporary solves to m0 instead of
 // allocating/factoring matrices whose answers cannot enter the iteration.
 int maximum_mode=-1;
 int mode_limit()const{return maximum_mode<0?modal->n[2]/2:maximum_mode;}
 int factor_group(int mode,int component)const{
  return 2*(mode<=mode_limit()?mode:0)+(component?1:0);
 }
 void factor(std::vector<ModalBlock>*vector_cache=nullptr){
  scalar_factorizations=vector_factorizations=0;
  if(modal){
   const int na=modal->n[0],nb=modal->n[1],half=mode_limit();
   blocks.resize(2*(half+1));
   const bool reuse=vector_cache&&!vector_cache->empty();
   if(reuse&&vector_cache->size()!=size_t(half+1))throw std::runtime_error("Invalid modal vector cache");
   if(reuse)for(const auto&B:*vector_cache)
    if(B.na!=na||B.nb!=nb||B.lu.size()!=size_t(na*na*nb)||B.transfer.size()!=size_t(na*na*nb)||
       B.lower.size()!=size_t(na*nb*(HISPID_CURVED_MODAL_PRECONDITIONER?3:1))||B.upper.size()!=size_t(na*nb*(HISPID_CURVED_MODAL_PRECONDITIONER?3:1))||B.permutation.size()!=size_t(na*nb))
     throw std::runtime_error("Incomplete modal vector cache");
   for(int mode=0;mode<=half;mode++)for(int v=0;v<2;v++){
    auto&B=blocks[2*mode+v];B.na=na;B.nb=nb;
    if(HISPID_FOCUS_TAU&&mode<=4)B.focus_endpoint=hispid::focus_weights(na,nb,mode);
    if(HISPID_EXACT_POLAR_TAU&&mode>=5)B.polar_endpoint=hispid::tau_weights(na,nb);
    if(v&&reuse){
     B=std::move((*vector_cache)[mode]);
     if(B.na!=na||B.nb!=nb)throw std::runtime_error("Modal vector cache grid mismatch");
     continue;
    }
    B.lu.assign(na*na*nb,0);B.transfer.assign(na*na*nb,0);
    B.banded_coupling=bool(HISPID_CURVED_MODAL_PRECONDITIONER);B.lower.assign(na*nb*(B.banded_coupling?3:1),0);B.upper.assign(na*nb*(B.banded_coupling?3:1),0);B.permutation.resize(na*nb);
    for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
     const int row=4*(i+na*(j+nb*mode))+v;
     for(size_t q=0;q<col[row].size();q++){
      const int local=col[row][q]/4-na*nb*mode,ci=local%na,cj=local/na;
      if(cj==j)B.lu[(j*na+i)*na+ci]+=val[row][q];
      else if(B.banded_coupling&&std::abs(ci-i)<=1&&cj==j-1)B.lower[3*(j*na+i)+ci-i+1]+=val[row][q];
      else if(B.banded_coupling&&std::abs(ci-i)<=1&&cj==j+1)B.upper[3*(j*na+i)+ci-i+1]+=val[row][q];
      else if(ci==i&&cj==j-1)B.lower[j*na+i]+=val[row][q];
      else if(ci==i&&cj==j+1)B.upper[j*na+i]+=val[row][q];
      else throw std::runtime_error("Modal FD stencil is not block tridiagonal");
     }
    }
    B.factor(2*mode+v);
    if(v)vector_factorizations++;else scalar_factorizations++;
   }
   if(vector_cache)vector_cache->clear();
   return;
  }
  int n=col.size();diag.resize(n);
  for(int i=0;i<n;i++){
   auto&cc=col[i];auto&vv=val[i];int di=std::lower_bound(cc.begin(),cc.end(),i)-cc.begin();diag[i]=di;
   for(int a=0;a<di;a++){
    int j=cc[a];double piv=val[j][diag[j]];if(std::abs(piv)<1e-30)throw std::runtime_error("ILU zero pivot");
    vv[a]/=piv;
    for(size_t q=diag[j]+1;q<col[j].size();q++){
     auto it=std::lower_bound(cc.begin()+a+1,cc.end(),col[j][q]);if(it!=cc.end()&&*it==col[j][q])vv[it-cc.begin()]-=vv[a]*val[j][q];
    }
   }
   if(!std::isfinite(vv[di])||std::abs(vv[di])<1e-30)throw std::runtime_error("ILU singular pivot");
  }
 }
 void retain_vectors(std::vector<ModalBlock>&cache){
  if(!modal)throw std::runtime_error("Vector reuse requires modal blocks");
  const int half=mode_limit();cache.resize(half+1);
  for(int mode=0;mode<=half;mode++)cache[mode]=std::move(blocks[2*mode+1]);
 }
 void project_rhs(const double*b,double*x)const{
   if(!modal)throw std::runtime_error("Modal projection requires modal factors");
   const int N=modal->n[2],stride=modal->n[0]*modal->n[1];
   for(int line=0;line<stride;line++)for(int v=0;v<4;v++){
    double sum=0,error=0;for(int k=0;k<N;k++){double y=b[4*(line+k*stride)+v]-error,t=sum+y;error=(t-sum)-y;sum=t;}
    const double mean=sum/N;
    for(int mode=0;mode<N;mode++){
     double value=mode==0?sum/std::sqrt(double(N)):0;
     if(mode>0){
      if(HISPID_COMPENSATED_MODAL_PROJECTION)value=puncture::compensated_projection(b+4*line+v,4*stride,modal->forward.data()+mode*N,N,mean);
      else for(int k=0;k<N;k++)value+=modal->forward[mode*N+k]*(b[4*(line+k*stride)+v]-mean);
     }
     const int row=4*(line+mode*stride)+v;x[row]=value/row_scale[row];
    }
   }
 }
 void solve(const double*b,double*x)const{
  int n=col.size();std::copy(b,b+n,x);
  if(modal){
   project_rhs(b,x);
   const int N=modal->n[2],half=N/2;
   for(int k=0;k<N;k++)for(int v=0;v<4;v++)
    blocks[factor_group(k<=half?k:k-half,v)].solve(x,k,v);
   return;
  }
  for(int i=0;i<n;i++)for(int k=0;k<diag[i];k++)x[i]-=val[i][k]*x[col[i][k]];
  for(int i=n-1;i>=0;i--){for(size_t k=diag[i]+1;k<col[i].size();k++)x[i]-=val[i][k]*x[col[i][k]];x[i]/=val[i][diag[i]];}
 }
 void solve(const std::vector<double>&b,std::vector<double>&x)const{x.resize(b.size());solve(b.data(),x.data());}
};
// Modal flat-Laplacian block preconditioner for the regular P unknowns.
// B_rm follows the analytic prolate Laplacian. Curved coefficients are
// represented by their azimuthal average; the vector block uses4/3 Delta.
struct AzimuthalAverage { double mu=0,potential=0; };
AzimuthalAverage azimuthal_average(const HiSpID_Data&s,int i,int j){
 AzimuthalAverage out;const int np=s.local.n[2];
 for(int phi=0;phi<np;phi++){
  int q=pindex(s,i,j,phi);const auto&g=s.geometry[q];const auto&u=s.basefields[q];double A[9];
  for(int x=0;x<9;x++){A[x]=g.M[x];for(int v=0;v<3;v++)for(int d=0;d<4;d++)A[x]+=g.L[x][v][d]*u[v+1][d];}
  const double psi=g.psi+u[0][0];out.mu+=(g.inv[0]+g.inv[4]+g.inv[8])/(3*np);
  out.potential+=g.g*(-g.R/8-5*std::pow(psi,4)*g.K*g.K/12-7*contraction(g,A,A)/(8*std::pow(psi,8)))/np;
 }
 return out;
}
Sparse preconditioner(HiSpID_Data&s,std::vector<ModalBlock>*vector_cache=nullptr,bool share_averages=true){
 if(s.quadrature){
#ifdef PUNCTURES_KOKKOS
  if(s.device)s.device->fields(s.values.data(),true);else
#endif
  fields(s,s.values.data(),s.basefields);
 }
 bool compact=bool(HISPID_SPECTRAL_RADIAL_PRECONDITIONER);
#ifdef PUNCTURES_KOKKOS
 puncture::TimedIf setup_timer(puncture::statistics().setup_seconds,bool(s.device));
 compact=compact||bool(s.device);if(compact)share_averages=true;
#endif
 Sparse mat;if(!compact){mat.col.resize(s.ntotal);mat.val.resize(s.ntotal);}mat.row_scale.resize(s.ntotal);mat.modal=&s.derivatives;
 const int na=s.local.n[0],nb=s.local.n[1],np=s.local.n[2],half=np/2;
 mat.maximum_mode=s.axisymmetric?1:half;
 const int factor_half=mat.mode_limit();
 const hispid::RadialPreconditioner radial(HISPID_SPECTRAL_RADIAL_PRECONDITIONER?na:0,hispid::AxisDerivatives::radial_stretch);
 const double ha=Pi/na,hb=Pi/nb;
 const auto endpoint=HISPID_AXIS_TAU?hispid::tau_weights(na,nb):std::vector<double>{};
 std::vector<AzimuthalAverage>averages;
 if(share_averages){
  averages.resize(na*nb);
#ifdef PUNCTURES_KOKKOS
  if(s.device){static_assert(sizeof(AzimuthalAverage)==2*sizeof(double));s.device->azimuthal_averages(reinterpret_cast<double*>(averages.data()));}else
#endif
  for(int j=0;j<nb;j++)for(int i=0;i<na;i++)averages[i+na*j]=azimuthal_average(s,i,j);
 }
 auto&curved=s.curved_symbols;
 if(HISPID_CURVED_MODAL_PRECONDITIONER&&curved.empty()){
  curved.resize(2*(half+1)*na*nb);
#ifdef PUNCTURES_KOKKOS
  if(s.device)s.device->curved_symbols(s.b,curved);else
#endif
  for(int mode=0;mode<=half;mode++)for(int group=0;group<2;group++)for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
   auto&avg=curved[(2*mode+group)*na*nb+i+na*j];
   for(int k=0;k<np;k++){
    const double A=s.derivatives.coordinate[0][i],B=s.derivatives.coordinate[1][j];
    double a=.5*(A+1),t=a*a,C=(1+t)/(1-t),eta=-2*B/(1+B*B),rho=s.b*std::sqrt((C*C-1)*(1-eta*eta)),phi=2*Pi*k/np;
    auto v=hispid::curved_symbol(s.geometry[pindex(s,i,j,k)],s.b,s.b*C*eta,rho*std::cos(phi),rho*std::sin(phi),mode,group);
    avg.zero+=v.zero/np;avg.mixed+=v.mixed/np;for(int d=0;d<2;d++){avg.first[d]+=v.first[d]/np;avg.second[d]+=v.second[d]/np;}
   }
  }
 }
 const bool reuse=compact&&vector_cache&&!vector_cache->empty();
 if(compact){mat.blocks.resize(2*(factor_half+1));if(reuse&&vector_cache->size()!=size_t(factor_half+1))throw std::runtime_error("invalid compact vector factors");
  for(int mode=0;mode<=factor_half;mode++)for(int v=0;v<2;v++){auto&B=mat.blocks[2*mode+v];
   if(v&&reuse){B=std::move((*vector_cache)[mode]);continue;}
   B.na=na;B.nb=nb;
    if(HISPID_FOCUS_TAU&&mode<=4)B.focus_endpoint=hispid::focus_weights(na,nb,mode);
    if(HISPID_EXACT_POLAR_TAU&&mode>=5)B.polar_endpoint=hispid::tau_weights(na,nb);B.lu.assign(size_t(na)*na*nb,0);B.transfer.assign(size_t(na)*na*nb,0);B.banded_coupling=bool(HISPID_CURVED_MODAL_PRECONDITIONER);B.lower.assign(na*nb*(B.banded_coupling?3:1),0);B.upper.assign(na*nb*(B.banded_coupling?3:1),0);B.permutation.resize(na*nb);
  }
 }
 for(int k=0;k<(compact?half+1:np);k++)for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
  const int mode=k<=half?k:k-half,r=hispid::AxisDerivatives::exponent(mode),p=pindex(s,i,j,k);
  const double a=.5*(s.derivatives.coordinate[0][i]+1),t=a*a,B=s.derivatives.coordinate[1][j],eta=-2*B/(1+B*B),s2=1-eta*eta,h2=4*t/((1-t)*(1-t)),D=s.b*s.b*(h2+s2);
  const double c=-2*(1-t)*std::pow(a*std::sqrt(s2),r),weight=row_weight(s,pindex(s,i,j,0));
  const auto average=share_averages?averages[i+na*j]:azimuthal_average(s,i,j);
  const double mu=average.mu,potential=average.potential;
  const double ctt=t*(1-t)*(1-t),ct=(1-t)*((r+1)*(1-t)-2*t),cee=s2,ce=-2*(r+1)*eta;
  const double c0=-(r+1)*(r+1-t)+(r*r-mode*mode)*(1/h2+1/s2);
  std::map<int,double>stencil;
  auto add=[&](int di,int dj,double v){int q=Index(0,i+di,j+dj,0,1,na,nb,1);stencil[q]+=v;};
  add(0,0,c0);
  const double al=ha*(i+.5),be=hb*(j+.5),sa=std::sin(al),sb=std::sin(be);
  const double sigma=.5*(1-std::cos(al));double dt,ddt;hispid::AxisDerivatives::radial_forward(sigma,dt,ddt);
  const double ta=.5*sa*dt,taa=.5*std::cos(al)*dt+.25*sa*sa*ddt;
  const double kappa=hispid::AxisDerivatives::angular_stretch,zeta=-std::cos(be),T=std::tanh(kappa),tanh=std::tanh(kappa*zeta),de=kappa*(1-tanh*tanh)/T,dde=-2*kappa*kappa*tanh*(1-tanh*tanh)/T,eb=de*sb,ebb=de*std::cos(be)+dde*sb*sb;
  const double ab=ct/ta-ctt*taa/(ta*ta*ta),bb=ce/eb-cee*ebb/(eb*eb*eb);
  double aa=ctt/(ta*ta),ba=cee/(eb*eb);
#if HISPID_MONOTONE_PRECONDITIONER
  // Minimal added diffusion makes both off-diagonal drift/diffusion weights
  // nonnegative. This changes only the approximate inverse, not the spectral
  // residual or its JVP. High regularity powers otherwise give cell Pe > 1.
  aa=std::max(aa,std::abs(ab)*ha/2);ba=std::max(ba,std::abs(bb)*hb/2);
#endif
  if(HISPID_SPECTRAL_RADIAL_PRECONDITIONER){
   // Off-diagonal differences preserve the constant radial mode exactly.
   double diagonal=0;for(int ii=0;ii<na;ii++)if(ii!=i){
    const double value=radial.entry(i,ii,r);
    add(ii-i,0,value);diagonal-=value;
   }add(0,0,diagonal);
  }else{
   add(1,0,aa/(ha*ha)+ab/(2*ha));add(-1,0,aa/(ha*ha)-ab/(2*ha));add(0,0,-2*aa/(ha*ha));
  }
  add(0,1,ba/(hb*hb)+bb/(2*hb));add(0,-1,ba/(hb*hb)-bb/(2*hb));add(0,0,-2*ba/(hb*hb));
  const bool tau_boundary=HISPID_AXIS_TAU&&mode>=5&&(i==0||j==0||j==nb-1);
  if(tau_boundary){
   stencil.clear();
   if(j==0||j==nb-1){
    int next=j==0?1:nb-2;double z0=-std::cos(hb*(j+.5)),z1=-std::cos(hb*(next+.5)),end=j==0?-1:1;
    add(0,0,(z1-end)/(z1-z0));add(0,next-j,(end-z0)/(z1-z0));
   }else{
    // The radial blocks are already dense: enforce the exact endpoint row
    // at no additional matrix storage, instead of a two-point approximation.
    for(int ii=0;ii<na;ii++)add(ii,0,endpoint[2*nb+ii]);
   }
  }
  for(int v=0;v<4;v++){
   const int row=4*p+v;
   std::map<int,double>curved_stencil;
   if(HISPID_CURVED_MODAL_PRECONDITIONER&&!tau_boundary){
    const auto&q=curved[(2*mode+(v?1:0))*na*nb+i+na*j];
    // Equilibrate the radial blocks geometrically; the opposite factor
    // in row_scale preserves the same physical approximate inverse.
    auto insert=[&](int di,int dj,double value){curved_stencil[Index(0,i+di,j+dj,0,1,na,nb,1)]+=D*value;};
    const double sa1=.5*sa,sa2=.5*std::cos(al),zb1=sb,zb2=std::cos(be);
    double ca=q.second[0]/(sa1*sa1),da=q.first[0]/sa1-q.second[0]*sa2/(sa1*sa1*sa1);
    double cb=q.second[1]/(zb1*zb1),db=q.first[1]/zb1-q.second[1]*zb2/(zb1*zb1*zb1);
#if HISPID_MONOTONE_PRECONDITIONER
    ca=std::max(ca,std::abs(da)*ha/2);cb=std::max(cb,std::abs(db)*hb/2);
#endif
    insert(0,0,q.zero+(v?0:potential));
    if(HISPID_SPECTRAL_RADIAL_PRECONDITIONER){
     double diagonal=0;for(int ii=0;ii<na;ii++)if(ii!=i){double z=4*q.second[0]*radial.second[i*na+ii]+2*q.first[0]*radial.first[i*na+ii];insert(ii-i,0,z);diagonal-=z;}insert(0,0,diagonal);
    }else{insert(1,0,ca/(ha*ha)+da/(2*ha));insert(-1,0,ca/(ha*ha)-da/(2*ha));insert(0,0,-2*ca/(ha*ha));}
    const double mixed=q.mixed/(sa1*zb1*4*ha*hb);
    insert(1,1,mixed);insert(-1,-1,mixed);insert(1,-1,-mixed);insert(-1,1,-mixed);
    insert(0,1,cb/(hb*hb)+db/(2*hb));insert(0,-1,cb/(hb*hb)-db/(2*hb));insert(0,0,-2*cb/(hb*hb));
   }
   const bool boundary=tau_boundary||(HISPID_FOCUS_TAU&&mode<=4&&i<hispid::focus_rank/2&&(j==0||j==nb-1));
   mat.row_scale[row]=boundary?1:(HISPID_CURVED_MODAL_PRECONDITIONER?weight*c/D:weight*mu*(v?4./3:1.)*c/D);
   const auto&entries=HISPID_CURVED_MODAL_PRECONDITIONER&&!tau_boundary?curved_stencil:stencil;
   for(auto entry:entries){int col=4*(entry.first+na*nb*k)+v;double value=entry.second;
    if(!HISPID_CURVED_MODAL_PRECONDITIONER&&col==row&&v==0&&!tau_boundary)value+=potential*D/mu;
    if(compact){
     if(k>factor_half||v>1||(v&&reuse))continue;
     auto&block=mat.blocks[2*k+v];int ci=entry.first%na,cj=entry.first/na;
     if(cj==j)block.lu[(j*na+i)*na+ci]+=value;
     else if(block.banded_coupling&&std::abs(ci-i)<=1&&cj==j-1)block.lower[3*(j*na+i)+ci-i+1]+=value;
     else if(block.banded_coupling&&std::abs(ci-i)<=1&&cj==j+1)block.upper[3*(j*na+i)+ci-i+1]+=value;
     else if(ci==i&&cj==j-1)block.lower[j*na+i]+=value;
     else if(ci==i&&cj==j+1)block.upper[j*na+i]+=value;
     else throw std::runtime_error("compact modal stencil is not block tridiagonal");continue;
    }
    if(value!=0||col==row){mat.col[row].push_back(col);mat.val[row].push_back(value);}
   }
  }
 }
 if(compact){
  for(int k=half+1;k<np;k++)for(int row=0;row<4*na*nb;row++)mat.row_scale[row+4*na*nb*k]=mat.row_scale[row+4*na*nb*(k-half)];
  auto factor_one=[&](int group){if(!(group%2)||!reuse)mat.blocks[group].factor(group);};
#ifdef PUNCTURES_KOKKOS
  if(s.device){
   std::mutex lock;std::string error;
   Kokkos::parallel_for("host modal factorizations",Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0,2*(factor_half+1)),[&](int group){try{factor_one(group);}catch(const std::exception&e){std::lock_guard<std::mutex>guard(lock);if(error.empty())error=e.what();}});
   Kokkos::DefaultHostExecutionSpace().fence();if(!error.empty())throw std::runtime_error(error);
  }else
#endif
  for(int group=0;group<2*(factor_half+1);group++)factor_one(group);
  mat.scalar_factorizations=factor_half+1;mat.vector_factorizations=reuse?0:factor_half+1;
  if(vector_cache)vector_cache->clear();
 }else mat.factor(vector_cache);return mat;
}
int relaxation_sweeps(const HiSpID_Data&s){return HISPID_DAMPED_MODAL_RELAXATION&&s.seed_family==HISPID_SEED_TRUMPET_R0_M?2:1;}
struct HiLinearContext {
 HiSpID_Data&data;const Sparse&preconditioner;
 std::vector<double>residual,action,correction,conditioned;
 explicit HiLinearContext(HiSpID_Data&s,const Sparse&M):data(s),preconditioner(M){
  bool host=true;
#ifdef PUNCTURES_KOKKOS
  host=!s.device;
#endif
  if(host&&(relaxation_sweeps(s)>1||(HISPID_LEFT_GMRES&&s.seed_family==HISPID_SEED_TRUMPET_R0_M))){residual.resize(s.ntotal);action.resize(s.ntotal);correction.resize(s.ntotal);if(HISPID_LEFT_GMRES)conditioned.resize(s.ntotal);}
 }
};
void hi_linear_progress(void*context,int iteration,double residual,double,double,double){
 if(iteration%(HISPID_CONSTRAINT_BLOCK_PRECONDITIONER?10:100))return;
 double rhs=*static_cast<double*>(context);
 std::fprintf(stderr,"HiSpID linear it=%d recurrence_relative=%.17g (true residual checked by controller)\n",iteration,rhs?residual/rhs:residual);std::fflush(stderr);
}
#ifdef PUNCTURES_KOKKOS
// HiSpID-only modal Woodbury update. Shared Bowen-York machinery is unchanged.
struct TauPolarDevice {
 int na=0,nb=0,np=0,groups=0;
 puncture::View response,inverse,delta,g,z,stream_rhs,stream_correction;
 puncture::Indices map;
 explicit TauPolarDevice(const Sparse&M){
  na=M.modal->n[0];nb=M.modal->n[1];np=M.modal->n[2];
  std::vector<int>indices(M.blocks.size(),-1),tasks(4*np);
  for(size_t i=0;i<M.blocks.size();i++)if(!M.blocks[i].polar_inverse.empty())indices[i]=groups++;
  if(!groups)return;
  const size_t rank=2*na,wide=size_t(na)*nb*rank,square=rank*rank;
  if(!HISPID_STREAM_POLAR_PRECONDITIONER)response=puncture::View(Kokkos::view_alloc(Kokkos::WithoutInitializing,"polar responses"),groups*wide);
  else{stream_rhs=puncture::View("implicit polar RHS",4ull*na*nb*np);stream_correction=puncture::View("implicit polar correction",4ull*na*nb*np);}
  inverse=puncture::View(Kokkos::view_alloc(Kokkos::WithoutInitializing,"polar inverses"),groups*square);
  delta=puncture::View(Kokkos::view_alloc(Kokkos::WithoutInitializing,"polar differences"),groups*2*nb);
  {puncture::Timed timer(puncture::statistics().transfer_seconds);
   for(size_t i=0;i<M.blocks.size();i++)if(indices[i]>=0){
    size_t k=indices[i];const auto&B=M.blocks[i];
    if(!HISPID_STREAM_POLAR_PRECONDITIONER)Kokkos::deep_copy(Kokkos::subview(response,std::make_pair(k*wide,(k+1)*wide)),puncture::Host(B.polar_response.data(),wide));
    Kokkos::deep_copy(Kokkos::subview(inverse,std::make_pair(k*square,(k+1)*square)),puncture::Host(B.polar_inverse.data(),square));
    Kokkos::deep_copy(Kokkos::subview(delta,std::make_pair(k*2*nb,(k+1)*2*nb)),puncture::Host(B.polar_delta.data(),2*nb));
   }
   puncture::statistics().host_to_device_bytes+=8*(response.extent(0)+inverse.extent(0)+delta.extent(0));
  }
  for(int k=0;k<np;k++)for(int v=0;v<4;v++)tasks[4*k+v]=indices[M.factor_group(k<=np/2?k:k-np/2,v)];
  map=puncture::Indices("polar task map",tasks.size());
  using H=Kokkos::View<const int*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  Kokkos::deep_copy(map,H(tasks.data(),tasks.size()));
  g=puncture::View("polar moments",4*np*rank);z=puncture::View("polar Schur solution",4*np*rank);
 }
 unsigned long long bytes()const{return 8*(response.extent(0)+inverse.extent(0)+delta.extent(0)+g.extent(0)+z.extent(0)+stream_rhs.extent(0)+stream_correction.extent(0))+4*map.extent(0);}
 void apply(puncture::View x,puncture::Modal*base=nullptr)const{
  if(!groups)return;
  const int a=na,b=nb,rank=2*a,tasks=4*np,N=a*b;
  auto W=response,S=inverse,V=delta,G=g,Z=z;auto mapping=map;
  Kokkos::parallel_for("polar boundary moments",Kokkos::RangePolicy<puncture::Exec>(0,tasks*rank),KOKKOS_LAMBDA(int q){
   int task=q/rank,k=q%rank,group=mapping(task);if(group<0)return;
   int side=k/a,i=k%a,mode=task/4,v=task%4;double sum=0;
   for(int j=0;j<b;j++)sum+=V((group*2+side)*b+j)*x(4*(i+a*j+N*mode)+v);
   G(q)=sum;
  });
  Kokkos::parallel_for("polar Schur solve",Kokkos::RangePolicy<puncture::Exec>(0,tasks*rank),KOKKOS_LAMBDA(int q){
   int task=q/rank,i=q%rank,group=mapping(task);if(group<0)return;double sum=0;
   for(int k=0;k<rank;k++)sum+=S((size_t(group)*rank+i)*rank+k)*G(task*rank+k);
   Z(q)=sum;
  });
  if(HISPID_STREAM_POLAR_PRECONDITIONER){
   if(!base)throw std::runtime_error("implicit polar correction requires its base inverse");
   auto rhs=stream_rhs,correction=stream_correction;
   Kokkos::parallel_for("implicit polar boundary injection",puncture::Range(0,tasks*N),KOKKOS_LAMBDA(int q){
    int task=q/N,row=q%N,i=row%a,j=row/a,group=mapping(task);double value=0;
    if(group>=0&&(j==0||j==b-1))value=Z(task*rank+(j? a:0)+i);
    rhs(4*(row+N*(task/4))+task%4)=value;
   });
   base->apply(rhs,correction,true);
   puncture::each(x,KOKKOS_LAMBDA(int i){x(i)-=correction(i);});return;
  }
  Kokkos::parallel_for("polar boundary correction",Kokkos::RangePolicy<puncture::Exec>(0,tasks*N),KOKKOS_LAMBDA(int q){
   int task=q/N,row=q%N,group=mapping(task);if(group<0)return;double sum=0;
   for(int k=0;k<rank;k++)sum+=W((size_t(group)*N+row)*rank+k)*Z(task*rank+k);
   x(4*(row+N*(task/4))+task%4)-=sum;
  });
 }
};
// Cartesian compatibility rows use the same exact boundary solve on host/device.
struct FocusDevice {
 int a=0,b=0,np=0,groups=0;puncture::View weight,response,inverse,moments,forward;puncture::Indices map;
 explicit FocusDevice(const Sparse&M):a(M.modal->n[0]),b(M.modal->n[1]),np(M.modal->n[2]){
  std::vector<int>indices(M.blocks.size(),-1),tasks(4*np);
  for(size_t k=0;k<M.blocks.size();k++)if(!M.blocks[k].focus_inverse.empty())indices[k]=groups++;
  if(!groups)return;const int N=a*b;
  std::vector<double>w,r,inv;
  for(const auto&B:M.blocks)if(!B.focus_inverse.empty()){
   w.insert(w.end(),B.focus_endpoint.begin(),B.focus_endpoint.end());r.insert(r.end(),B.focus_response.begin(),B.focus_response.end());inv.insert(inv.end(),B.focus_inverse.begin(),B.focus_inverse.end());
  }
  weight=puncture::upload(w,"focus border weights");response=puncture::upload(r,"focus responses");inverse=puncture::upload(inv,"focus inverses");forward=puncture::upload(M.modal->forward,"focus Fourier projection");
  for(int k=0;k<np;k++)for(int v=0;v<4;v++)tasks[4*k+v]=indices[M.factor_group(k<=np/2?k:k-np/2,v)];
  map=puncture::Indices("focus task map",tasks.size());using H=Kokkos::View<const int*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;Kokkos::deep_copy(map,H(tasks.data(),tasks.size()));
  moments=puncture::View("focus boundary defects",4*hispid::focus_rank*np);
 }
 unsigned long long bytes()const{return 8*(weight.extent(0)+response.extent(0)+inverse.extent(0)+moments.extent(0)+forward.extent(0))+4*map.extent(0);}
 void apply(puncture::View rhs,puncture::View x,bool projected=false)const{
  if(!groups)return;const int N=a*b,P=np,A=a,B=b;
  auto W=weight,R=response,S=inverse,G=moments,F=forward;auto mapping=map;
  using Team=Kokkos::TeamPolicy<puncture::Exec>;using Member=Team::member_type;
  Kokkos::parallel_for("focus inverse moments",Team(4*hispid::focus_rank*P,puncture::reduction_team_size),KOKKOS_LAMBDA(const Member&team){
   int q=team.league_rank(),task=q/hispid::focus_rank,side=q%hispid::focus_rank,group=mapping(task),mode=task/4,v=task%4;if(group<0)return;
   puncture::ProjectionAccumulator value;Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,N),[&](int p,puncture::ProjectionAccumulator&sum){sum.product(W((hispid::focus_rank*group+side)*N+p),x(4*(p+N*mode)+v));},Kokkos::Sum<puncture::ProjectionAccumulator>(value));
   Kokkos::single(Kokkos::PerTeam(team),[&](){int row=hispid::focus_row(side,A,B);double source;
    if(projected)source=rhs(4*(row+N*mode)+v);
    else{double sum=0,error=0;for(int k=0;k<P;k++){double y=rhs(4*(row+N*k)+v)-error,t=sum+y;error=(t-sum)-y;sum=t;}
     source=mode?0:sum/Kokkos::sqrt(double(P));
     if(mode){if(HISPID_COMPENSATED_MODAL_PROJECTION)source=puncture::compensated_projection(rhs.data()+4*row+v,4*N,F.data()+mode*P,P,sum/P);else for(int k=0;k<P;k++)source+=F(mode*P+k)*(rhs(4*(row+N*k)+v)-sum/P);}}
    G(q)=source-(value.high+value.low);
   });
  });
  Kokkos::parallel_for("focus inverse correction",puncture::Range(0,4*P*N),KOKKOS_LAMBDA(int q){
   int task=q/N,p=q%N,group=mapping(task);if(group<0)return;double change=0;
   for(int i=0;i<hispid::focus_rank;i++)for(int j=0;j<hispid::focus_rank;j++)change+=R((group*N+p)*hispid::focus_rank+i)*S(hispid::focus_rank*hispid::focus_rank*group+hispid::focus_rank*i+j)*G(hispid::focus_rank*task+j);
   x(4*(p+N*(task/4))+task%4)+=change;
  });
 }
};
puncture::Modal device_preconditioner(const Sparse&M){
 const int a=M.modal->n[0],b=M.modal->n[1],p=M.modal->n[2],groups=M.blocks.size();
 puncture::Modal d(a,b,p,4,groups,true,true);d.lower_bandwidth=HISPID_CURVED_MODAL_PRECONDITIONER?1:0;
 d.compensated_projection=bool(HISPID_COMPENSATED_MODAL_PROJECTION);d.skip_zero_rhs=bool(HISPID_CONSTRAINT_BLOCK_PRECONDITIONER);
 std::vector<int>rows(4*p);
 if constexpr(Kokkos::SpaceAccessibility<puncture::Exec,Kokkos::HostSpace>::accessible){
  std::vector<puncture::ModalPointers>ptrs;for(const auto&B:M.blocks)ptrs.push_back({B.lu.data(),B.transfer.data(),B.lower.data(),B.permutation.data()});
  d.native_blocks=decltype(d.native_blocks)("borrowed modal factors",groups);
  using PH=Kokkos::View<const puncture::ModalPointers*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  Kokkos::deep_copy(d.native_blocks,PH(ptrs.data(),groups));
  d.row_scale=puncture::View(const_cast<double*>(M.row_scale.data()),M.row_scale.size());d.forward=puncture::View(const_cast<double*>(M.modal->forward.data()),M.modal->forward.size());
 }else{
  size_t matrix=size_t(a)*a*b,lower=size_t(a)*b*(HISPID_CURVED_MODAL_PRECONDITIONER?3:1);
  d.lu=puncture::View(Kokkos::view_alloc(Kokkos::WithoutInitializing,"modal LU"),groups*matrix);d.transfer=puncture::View(Kokkos::view_alloc(Kokkos::WithoutInitializing,"modal transfer"),groups*matrix);d.lower=puncture::View(Kokkos::view_alloc(Kokkos::WithoutInitializing,"modal lower"),groups*lower);
  std::vector<int>permutation;for(int g=0;g<groups;g++){const auto&B=M.blocks[g];
   {puncture::Timed t(puncture::statistics().transfer_seconds);
    Kokkos::deep_copy(Kokkos::subview(d.lu,std::make_pair(g*matrix,(g+1)*matrix)),puncture::Host(B.lu.data(),matrix));
    Kokkos::deep_copy(Kokkos::subview(d.transfer,std::make_pair(g*matrix,(g+1)*matrix)),puncture::Host(B.transfer.data(),matrix));
    Kokkos::deep_copy(Kokkos::subview(d.lower,std::make_pair(g*lower,(g+1)*lower)),puncture::Host(B.lower.data(),lower));puncture::statistics().host_to_device_bytes+=(2*matrix+lower)*8;
   }for(auto x:B.permutation)permutation.push_back(int(x));
  }
  d.permutation=puncture::Indices("modal pivots",permutation.size());using H=Kokkos::View<const int*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;Kokkos::deep_copy(d.permutation,H(permutation.data(),permutation.size()));
  d.row_scale=puncture::upload(M.row_scale,"modal scaling");d.forward=puncture::upload(M.modal->forward,"modal Fourier");
 }
 for(int k=0;k<p;k++)for(int v=0;v<4;v++)rows[4*k+v]=M.factor_group(k<=p/2?k:k-p/2,v);
 d.block=puncture::Indices("modal groups",rows.size());
 using H=Kokkos::View<const int*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
 Kokkos::deep_copy(d.block,H(rows.data(),rows.size()));d.prepare();return d;
}
#endif
int hi_linear_action(void*context,const double*input,double*output){
 auto&c=*static_cast<HiLinearContext*>(context);
 try{jvp(c.data,input,output);project_solve_sector(c.data,output,false);return 0;}catch(const std::exception&e){hispid::last_error=e.what();return -1;}
}
int hi_linear_precondition(void*context,const double*input,double*output){
 auto&c=*static_cast<HiLinearContext*>(context);const int n=c.data.ntotal;
 try{
  puncture::relaxed_inverse(input,output,c.residual.data(),c.action.data(),c.correction.data(),
   [&](auto in,auto out){c.data.preconditioner_applications++;c.preconditioner.solve(in,out);project_solve_sector(c.data,out,true);},
   [&](auto in,auto out){jvp(c.data,in,out);project_solve_sector(c.data,out,false);},
   [n](auto out,auto in){std::copy_n(in,n,out);},
   [n](auto out,double a){for(int i=0;i<n;i++)out[i]*=a;},
   [n](auto out,double a,auto in){for(int i=0;i<n;i++)out[i]+=a*in[i];},relaxation_sweeps(c.data),.4);
  return 0;
 }catch(const std::exception&e){hispid::last_error=e.what();return -1;}
}
#ifdef PUNCTURES_KOKKOS
// Fine-grid upper-triangular constraint-block inverse. The momentum block
// is independent of the scalar correction wherever grad K=0. Inner Krylov
// vectors store only the active component(s), without changing either block
// action or the original coupled stopping norm.
struct ConstraintBlockInverse {
 HiSpID_Data&data;puncture::Action inverse;
 puncture::View rhs,delta,action,input_full,output_full,packed_rhs,packed_solution;
 int applications=0,scalar_iterations=0,vector_iterations=0,scalar_incomplete=0,vector_incomplete=0;
 double last_scalar_relative=0,last_vector_relative=0;bool trace=false;
 // A single short inner cycle can return the same deficient vector space
 // on every outer application. Retain three previous-cycle corrections in
 // the shared LGMRES controller; acceptance still uses the original block
 // residual. The outer solve remains flexible for incomplete inner solves.
 static constexpr int limit=1024,restart=256,augmentation=3,workspace_bytes_per_point=208;
 ConstraintBlockInverse(HiSpID_Data&s,puncture::Action M):data(s),inverse(std::move(M)),
  rhs("constraint block RHS",s.ntotal),delta("constraint block correction",s.ntotal),action("constraint block action",s.ntotal),
  input_full("block full input",s.ntotal),output_full("block full output",s.ntotal),
  packed_rhs("packed block RHS",3ull*s.npt),packed_solution("packed block solution",3ull*s.npt){
  const char*setting=std::getenv("HISPID_TRACE_NEWTON");trace=setting&&std::strcmp(setting,"1")==0;
 }
 static void pack(puncture::View full,puncture::View packed,bool scalar){
  const int components=scalar?1:3;
  puncture::each(packed,KOKKOS_LAMBDA(int i){packed(i)=full(4*(i/components)+(scalar?0:i%components+1));});
 }
 static void unpack(puncture::View packed,puncture::View full,bool scalar){
  puncture::each(full,KOKKOS_LAMBDA(int i){int v=i%4,p=i/4;full(i)=scalar?(v==0?packed(p):0.):(v?packed(3*p+v-1):0.);});
 }
 struct Progress {double norm;bool scalar;int application;};
 static void progress(void*context,int iteration,double residual,double,double,double){
  if(iteration==0||iteration%restart)return;
  const auto&p=*static_cast<Progress*>(context);
  std::fprintf(stderr,"HiSpID inner block=%s application=%d it=%d recurrence_relative=%.17g (original true residual required)\n",
   p.scalar?"scalar":"vector",p.application,iteration,residual/p.norm);std::fflush(stderr);
 }
 void solve_block(puncture::View b,puncture::View x,bool scalar){
  const int n=(scalar?1:3)*data.npt;
  puncture::View bpack(packed_rhs.data(),n),xpack(packed_solution.data(),n);pack(b,bpack,scalar);
  const double norm=puncture::length(bpack);if(norm==0){puncture::zero(x);return;}
  PK_Options options{PK_LGMRES,limit,std::min(restart,n),1,0,0,.1*norm,nullptr,nullptr};PK_Result result{};Progress monitor{norm,scalar,applications};
  const int status=puncture::solve(bpack,xpack,options,
   [&](auto in,auto out){unpack(in,input_full,scalar);device_apply(data,input_full,output_full);pack(output_full,out,scalar);},
   [&](auto in,auto out){unpack(in,input_full,scalar);inverse(input_full,output_full);pack(output_full,out,scalar);},result,trace?progress:nullptr,&monitor,true);
  data.jvp_applications+=result.operator_calls;data.preconditioner_applications+=result.preconditioner_calls;
  (scalar?scalar_iterations:vector_iterations)+=result.iterations;
  (scalar?last_scalar_relative:last_vector_relative)=result.relative_residual;
  if(status!=PK_SUCCESS)(scalar?scalar_incomplete:vector_incomplete)++;
  if(status==PK_CALLBACK||status==PK_NONFINITE||status==PK_ALLOCATION||status==PK_INVALID||!puncture::finite(xpack)){
   hispid::last_error=std::string("constraint block inverse: ")+PK_status_string(status);throw std::runtime_error(hispid::last_error);
  }
  unpack(xpack,x,scalar);
 }
 void apply(puncture::View b,puncture::View x){
  applications++;solve_block(b,x,false);device_apply(data,x,action);data.jvp_applications++;
  auto r=rhs,a=action;puncture::each(r,KOKKOS_LAMBDA(int i){r(i)=i%4==0?b(i)-a(i):0.;});
  solve_block(rhs,delta,true);puncture::axpy(x,1.,delta);
  if(trace&&(applications==1||applications%10==0))report();
 }
 unsigned long long bytes()const{return 8*(rhs.extent(0)+delta.extent(0)+action.extent(0)+input_full.extent(0)+output_full.extent(0)+packed_rhs.extent(0)+packed_solution.extent(0));}
 void report()const{
  std::fprintf(stderr,"HiSpID constraint blocks applications=%d vector_iterations=%d scalar_iterations=%d vector_incomplete=%d scalar_incomplete=%d last_vector_relative=%.17g last_scalar_relative=%.17g inner_target=0.1 inner_limit=%d restart=%d method=lgmres\n",applications,vector_iterations,scalar_iterations,vector_incomplete,scalar_incomplete,last_vector_relative,last_scalar_relative,limit,restart);std::fflush(stderr);
 }
};
// Diagnostic of the frozen first VECTOR-inner problem, not a solve.
// Split modal unknowns, but measure actions in original momentum rows.
void audit_modal_response(HiSpID_Data&s,puncture::View input,const puncture::Action&inverse,const std::vector<double>&scales,puncture::View before_focus){
 if(s.quadrature)throw std::runtime_error("modal response diagnostic requires collocation execution");
 const int na=s.local.n[0],nb=s.local.n[1],np=s.local.n[2],N=na*nb;
 puncture::View b("modal audit vector RHS",s.ntotal),z("modal audit inverse",s.ntotal),part("modal audit split",s.ntotal),out("modal audit action",s.ntotal);
 puncture::each(b,KOKKOS_LAMBDA(int q){b(q)=q%4?input(q):0.;});
 const double norm=puncture::length(b);if(!(norm>0))throw std::runtime_error("empty modal audit vector RHS");
 puncture::scale_copy(b,1/norm,b);inverse(b,z);s.preconditioner_applications++;
 std::vector<double>rhs(s.ntotal),action(s.ntotal),combined(s.ntotal);std::vector<long double>sum(s.ntotal,0.);
 puncture::download(b,rhs.data());
 if(const char*prefix=std::getenv("HISPID_MODAL_RESPONSE_DUMP")){
  std::vector<double>correction(s.ntotal);puncture::download(z,correction.data());
  auto write=[&](const char*suffix,const std::vector<double>&values){
   std::string path=std::string(prefix)+suffix;FILE*file=std::fopen(path.c_str(),"wbx");
   if(!file)throw std::runtime_error("cannot create exclusive modal-response dump: "+path);
   const size_t count=std::fwrite(values.data(),sizeof(double),values.size(),file);const int closed=std::fclose(file);
   if(count!=values.size()||closed)throw std::runtime_error("incomplete modal-response dump: "+path);
  };
  write(".rhs.f64",rhs);write(".correction.f64",correction);write(".scale.f64",scales);
  if(before_focus.extent(0)){puncture::download(before_focus,correction.data());write(".before-focus.f64",correction);}
 }
 std::fprintf(stderr,"HiSpID modal response original_vector_rhs_l2=%.17g normalization=unit_vector_rhs groups=0,1,2,3,4,ge5\n",norm);
 for(int group=0;group<6;group++){
  puncture::each(part,KOKKOS_LAMBDA(int q){int slot=(q/4)/N,mode=slot<=np/2?slot:slot-np/2;part(q)=q%4&&(group<5?mode==group:mode>=5)?z(q):0.;});
  device_apply(s,part,out);s.jvp_applications++;auto f=s.device->op.fields;double value2=0,gradient2=0;
  Kokkos::parallel_reduce("modal audit physical value",puncture::Range(0,s.npt),KOKKOS_LAMBDA(int p,double&v){for(int c=1;c<4;c++)v+=f(40*p+10*c)*f(40*p+10*c);},value2);
  Kokkos::parallel_reduce("modal audit physical gradient",puncture::Range(0,s.npt),KOKKOS_LAMBDA(int p,double&v){for(int c=1;c<4;c++)for(int d=1;d<4;d++)v+=f(40*p+10*c+d)*f(40*p+10*c+d);},gradient2);
  puncture::download(out,action.data());for(int q=0;q<s.ntotal;q++)if(q%4)sum[q]+=action[q];
  long double input2=0,ordinary2=0,tau2=0,focus2=0;double maximum=0;int mi=-1,mj=-1,mm=-1,mv=-1;
  for(int slot=0;slot<np;slot++){
   int mode=slot<=np/2?slot:slot-np/2;
   for(int j=0;j<nb;j++)for(int i=0;i<na;i++)for(int v=1;v<4;v++){
    long double bv=0,av=0;for(int k=0;k<np;k++){int q=4*(i+na*j+N*k)+v;double w=s.derivatives.forward[slot*np+k];bv+=(long double)w*rhs[q];av+=(long double)w*action[q];}
    if(group<5?mode==group:mode>=5)input2+=bv*bv;
    bool focus_row=false;if(HISPID_FOCUS_TAU&&mode<=4)for(int side=0;side<hispid::focus_rank;side++)if(i+na*j==hispid::focus_row(side,na,nb))focus_row=true;
    const bool tau=HISPID_AXIS_TAU&&mode>=5&&(i==0||j==0||j==nb-1);
    if(focus_row)focus2+=av*av;else if(tau)tau2+=av*av;else{
     ordinary2+=av*av;if(std::abs(av)>maximum){maximum=std::abs(av);mi=i;mj=j;mm=mode;mv=v;}
    }
   }
  }
  double radius=std::numeric_limits<double>::quiet_NaN();
  if(mi>=0){const double a=.5*(s.derivatives.coordinate[0][mi]+1),t=a*a,C=(1+t)/(1-t),B=s.derivatives.coordinate[1][mj],eta=-2*B/(1+B*B),rho2=s.b*s.b*(C*C-1)*(1-eta*eta);
   int hole=s.local.hole[0].mass>0?0:1;double axial=s.b*C*eta-s.local.hole[hole].center[0];radius=std::sqrt(axial*axial+rho2);
  }
  std::fprintf(stderr,"HiSpID modal response group=%d rhs_group_l2=%.17g action_l2=%.17g pde_l2=%.17g focus_l2=%.17g tau_l2=%.17g physical_value_collocation_l2=%.17g physical_gradient_collocation_l2=%.17g pde_max=%.17g maximum_mode=%d maximum_component=%d maximum_i=%d maximum_j=%d maximum_radius=%.17g\n",group,double(std::sqrt(input2)),double(std::sqrt(ordinary2+focus2+tau2)),double(std::sqrt(ordinary2)),double(std::sqrt(focus2)),double(std::sqrt(tau2)),std::sqrt(value2),std::sqrt(gradient2),maximum,mm,mv,mi,mj,radius);std::fflush(stderr);
 }
 device_apply(s,z,out);s.jvp_applications++;puncture::download(out,combined.data());long double total=0,error=0;
 for(int q=0;q<s.ntotal;q++)if(q%4){total+=(long double)combined[q]*combined[q];long double e=sum[q]-combined[q];error+=e*e;}
 std::fprintf(stderr,"HiSpID modal response sum_action_l2=%.17g superposition_difference_l2=%.17g superposition_relative=%.17g jvp_count=7 krylov_iterations=0\n",double(std::sqrt(total)),double(std::sqrt(error)),double(std::sqrt(error/total)));std::fflush(stderr);
 if(before_focus.extent(0)){
  std::vector<double>base(s.ntotal),border(s.ntotal);auto prior=before_focus;
  for(int component=0;component<2;component++){
   puncture::each(part,KOKKOS_LAMBDA(int q){int slot=(q/4)/N,mode=slot<=np/2?slot:slot-np/2;part(q)=q%4&&mode==4?(component?z(q)-prior(q):prior(q)):0.;});
   device_apply(s,part,out);s.jvp_applications++;puncture::download(out,component?border.data():base.data());
  }
  long double aa=0,bb=0,ab=0,cc=0,af=0,bf=0,at=0,bt=0;
  for(int slot=0;slot<np;slot++)for(int j=0;j<nb;j++)for(int i=0;i<na;i++)for(int v=1;v<4;v++){
   int mode=slot<=np/2?slot:slot-np/2;long double av=0,bv=0;
   for(int k=0;k<np;k++){int q=4*(i+na*j+N*k)+v;double w=s.derivatives.forward[slot*np+k];av+=(long double)w*base[q];bv+=(long double)w*border[q];}
   bool boundary=false;if(HISPID_FOCUS_TAU&&mode<=4)for(int r=0;r<hispid::focus_rank;r++)if(i+na*j==hispid::focus_row(r,na,nb))boundary=true;
   bool tau=HISPID_AXIS_TAU&&mode>=5&&(i==0||j==0||j==nb-1);
   if(boundary){af+=av*av;bf+=bv*bv;}else if(tau){at+=av*av;bt+=bv*bv;}
   else{aa+=av*av;bb+=bv*bv;ab+=av*bv;cc+=(av+bv)*(av+bv);}
  }
  std::fprintf(stderr,"HiSpID focus response mode=4 base_ordinary_l2=%.17g border_ordinary_l2=%.17g ordinary_dot=%.17g combined_ordinary_l2=%.17g base_focus_rows_l2=%.17g border_focus_rows_l2=%.17g base_tau_l2=%.17g border_tau_l2=%.17g extra_jvp_count=2\n",double(std::sqrt(aa)),double(std::sqrt(bb)),double(ab),double(std::sqrt(cc)),double(std::sqrt(af)),double(std::sqrt(bf)),double(std::sqrt(at)),double(std::sqrt(bt)));std::fflush(stderr);
 }

}
bool constraint_blocks_enabled(const HiSpID_Data&s,int method){
 return HISPID_CONSTRAINT_BLOCK_PRECONDITIONER&&s.device&&!s.quadrature&&method!=PK_BICGSTAB&&
  !(HISPID_LEFT_GMRES&&s.seed_family==HISPID_SEED_TRUMPET_R0_M&&method==PK_GMRES);
}
#endif
bool linear_solve(HiSpID_Data&s,const Sparse&M,const std::vector<double>&rhs,
                  std::vector<double>&x,double rtol,int method,bool eager=false){
 x.assign(s.ntotal,0);HiLinearContext context{s,M};
 const char*trace_setting=std::getenv("HISPID_TRACE_NEWTON");const bool trace_linear=trace_setting&&std::strcmp(trace_setting,"1")==0;
 double rhs_norm=norm2v(rhs);PK_Monitor monitor=trace_linear?hi_linear_progress:nullptr;
 const char*action_probe_setting=std::getenv("HISPID_PROBE_LINEAR_ACTION");
 const bool action_probe=action_probe_setting&&std::strcmp(action_probe_setting,"1")==0;
 PK_Options options{method,s.local.max_krylov,s.local.krylov_restart,1,0,int(eager),rtol*norm2v(rhs),nullptr,nullptr};
 PK_Result result{};int status;
 const bool left=HISPID_LEFT_GMRES&&s.seed_family==HISPID_SEED_TRUMPET_R0_M&&method==PK_GMRES;
 auto progress=[&](int iteration,double relative){if(trace_linear){std::fprintf(stderr,"HiSpID left cycle it=%d original_relative=%.17g target=%.17g\n",iteration,relative,rtol);std::fflush(stderr);}};
#ifdef PUNCTURES_KOKKOS
 if(s.device){
  auto&Mstats=puncture::statistics();
  double begin=puncture::seconds();auto d=device_preconditioner(M);TauPolarDevice polar(M);FocusDevice focus(M);puncture::Exec().fence();Mstats.setup_seconds+=puncture::seconds()-begin;
  auto b=puncture::upload(rhs,"linear RHS");puncture::View solution("linear solution",s.ntotal);
  const char*modal_audit=std::getenv("HISPID_PROBE_MODAL_RESPONSE");
  if(modal_audit&&std::strcmp(modal_audit,"1")==0){
   puncture::View before_focus;
   const char*focus_audit=std::getenv("HISPID_PROBE_FOCUS_RESPONSE");
   if(focus_audit&&std::strcmp(focus_audit,"1")==0){
    if(!HISPID_FOCUS_TAU||s.local.n[2]/2<4)throw std::runtime_error("focus response audit requires focus rows and Fourier mode four");
    before_focus=puncture::View("modal audit before focus",s.ntotal);
    if(s.local.n[2]/2>=4){const auto&B=M.blocks[M.factor_group(4,1)];const int rank=hispid::focus_rank,N=s.local.n[0]*s.local.n[1];
     for(int i=0;i<rank;i++)for(int j=0;j<rank;j++){long double value=0;for(int p=0;p<N;p++)value+=(long double)B.focus_endpoint[i*N+p]*B.focus_response[rank*p+j];
      std::fprintf(stderr,"HiSpID focus Schur mode=4 row=%d col=%d value=%.17g inverse=%.17g\n",i,j,double(value),B.focus_inverse[rank*i+j]);}
    }
   }
   audit_modal_response(s,b,[&](auto in,auto out){d.apply(in,out);polar.apply(out,&d);if(before_focus.extent(0)){puncture::copy(before_focus,out);project_solve_sector(s,before_focus,true);}focus.apply(in,out);project_solve_sector(s,out,true);},M.row_scale,before_focus);
   s.last_gmres_relative=1.;hispid::last_error="modal response diagnostic only; no Krylov solve attempted";return false;
  }
  std::unique_ptr<ConstraintBlockInverse>blocks;
  if(constraint_blocks_enabled(s,method))blocks=std::make_unique<ConstraintBlockInverse>(s,[&](auto in,auto out){d.apply(in,out);polar.apply(out,&d);focus.apply(in,out);project_solve_sector(s,out,true);});
  puncture::View relaxation_residual,relaxation_action,relaxation_correction,left_conditioned;
  if(relaxation_sweeps(s)>1||left){
   relaxation_residual=puncture::View("relaxation residual",s.ntotal);relaxation_action=puncture::View("relaxation action",s.ntotal);relaxation_correction=puncture::View("relaxation correction",s.ntotal);
   if(left)left_conditioned=puncture::View("left conditioned RHS",s.ntotal);
  }
  if(action_probe){
   std::vector<double>reference,device(s.ntotal);puncture::View probe_rhs;M.solve(rhs,reference);
   project_solve_sector(s,reference.data(),true);
   for(int matched=0;matched<2;matched++){
   if(matched){
    std::vector<double>projected(s.ntotal);M.project_rhs(rhs.data(),projected.data());
    probe_rhs=puncture::upload(projected,"matched modal RHS");d.apply(probe_rhs,solution,true);
   }else d.apply(b,solution);
   polar.apply(solution,&d);focus.apply(matched?probe_rhs:b,solution,bool(matched));project_solve_sector(s,solution,true);
   puncture::download(solution,device.data());
   const int np=s.local.n[2],stride=s.local.n[0]*s.local.n[1];
   for(int mode=0;mode<=np/2;mode++)for(int v=0;v<4;v++){
    long double norm=0,error=0;double maximum=0;
    for(int slot=0;slot<np;slot++)if((slot<=np/2?slot:slot-np/2)==mode)
     for(int i=0;i<stride;i++){
      int q=4*(i+stride*slot)+v;long double delta=(long double)device[q]-reference[q];
      norm+=(long double)reference[q]*reference[q];error+=delta*delta;
      maximum=std::max(maximum,double(std::abs(delta)));
     }
    std::fprintf(stderr,"HiSpID %s inverse mode=%d component=%d reference_l2=%.17g difference_l2=%.17g difference_linf=%.17g\n",
      matched?"matched projection":"device",mode,v,double(std::sqrt(norm)),double(std::sqrt(error)),maximum);
   }
   }
   std::fflush(stderr);Kokkos::deep_copy(solution,0.);
  }
  Mstats.resident_bytes=std::max(Mstats.resident_bytes,(blocks?blocks->bytes():0)+s.device->bytes()+(s.quadrature?s.quadrature->device->bytes()+s.device_transfer->bytes():0)+polar.bytes()+focus.bytes()+8*(d.lu.extent(0)+d.inverse.extent(0)+d.transfer.extent(0)+d.lower.extent(0)+d.row_scale.extent(0)+d.forward.extent(0)+d.workspace.extent(0)+d.column.extent(0)+b.extent(0)+solution.extent(0)+relaxation_residual.extent(0)+relaxation_action.extent(0)+relaxation_correction.extent(0)+left_conditioned.extent(0)));
  if(left){
   status=puncture::left_gmres(b,solution,relaxation_residual,left_conditioned,relaxation_correction,relaxation_action,
    options,rtol,[](auto v){return puncture::length(v);},
    [](puncture::View out,puncture::View rhs,puncture::View ax){puncture::each(out,KOKKOS_LAMBDA(int i){out(i)=rhs(i)-ax(i);});},
    [](auto out,auto in){puncture::axpy(out,1.,in);},
    [&](auto in,auto out){device_apply(s,in,out);project_solve_sector(s,out,false);s.jvp_applications++;},
    [&](auto in,auto out){d.apply(in,out);polar.apply(out,&d);focus.apply(in,out);project_solve_sector(s,out,true);s.preconditioner_applications++;},
    [&](auto rhs,auto delta,const PK_Options&o,auto action,PK_Result&part,double scale){
     return puncture::solve(rhs,delta,o,action,[](auto in,auto out){puncture::copy(out,in);},part,monitor,&scale);
    },progress,result);
  }else
  status=puncture::solve(b,solution,options,[&](auto in,auto out){device_apply(s,in,out);project_solve_sector(s,out,false);},[&](auto in,auto out){
    if(blocks)blocks->apply(in,out);else
    puncture::relaxed_inverse(in,out,relaxation_residual,relaxation_action,relaxation_correction,
     [&](auto input,auto output){d.apply(input,output);polar.apply(output,&d);focus.apply(input,output);project_solve_sector(s,output,true);},
     [&](auto input,auto output){device_apply(s,input,output);project_solve_sector(s,output,false);},
     [](auto output,auto input){puncture::copy(output,input);},
     [](auto output,double a){puncture::scale_copy(output,a,output);},
     [](auto output,double a,auto input){puncture::axpy(output,a,input);},relaxation_sweeps(s),.4);
   },result,monitor,&rhs_norm);
  if(blocks&&trace_linear)blocks->report();
  puncture::download(solution,x.data());if(!left){s.jvp_applications+=result.operator_calls+(blocks?0:(relaxation_sweeps(s)-1)*result.preconditioner_calls);s.preconditioner_applications+=(blocks?0:relaxation_sweeps(s))*result.preconditioner_calls;}
 }else
#endif
 {
  if(left){
   status=puncture::left_gmres(const_cast<double*>(rhs.data()),x.data(),context.residual.data(),context.conditioned.data(),context.correction.data(),context.action.data(),
    options,rtol,[&](auto v){double q=0;for(int i=0;i<s.ntotal;i++)q+=v[i]*v[i];return std::sqrt(q);},
    [&](auto out,auto rhs,auto ax){for(int i=0;i<s.ntotal;i++)out[i]=rhs[i]-ax[i];},
    [&](auto out,auto in){for(int i=0;i<s.ntotal;i++)out[i]+=in[i];},
    [&](auto in,auto out){jvp(s,in,out);project_solve_sector(s,out,false);},
    [&](auto in,auto out){M.solve(in,out);project_solve_sector(s,out,true);s.preconditioner_applications++;},
    [&](auto rhs,auto delta,const PK_Options&o,auto action,PK_Result&part,double scale){
     using Action=decltype(action);
     struct Callbacks {
      Action&action;int n;
      static int apply(void*c,const double*in,double*out){try{static_cast<Callbacks*>(c)->action(const_cast<double*>(in),out);return 0;}catch(const std::exception&e){hispid::last_error=e.what();return -1;}}
      static int identity(void*c,const double*in,double*out){std::copy_n(in,static_cast<Callbacks*>(c)->n,out);return 0;}
     } callbacks{action,s.ntotal};
     return PK_solve(s.ntotal,rhs,delta,&o,Callbacks::apply,Callbacks::identity,&callbacks,monitor,&scale,&part);
    },progress,result);
  }else status=PK_solve(s.ntotal,rhs.data(),x.data(),&options,hi_linear_action,hi_linear_precondition,&context,monitor,&rhs_norm,&result);
 }
 s.diag.krylov_iterations+=result.iterations;s.last_gmres_relative=result.relative_residual;
 if(action_probe){
  std::vector<double>action(s.ntotal);jvp(s,x.data(),action.data());project_solve_sector(s,action.data(),false);
  const int na=s.local.n[0],nb=s.local.n[1],np=s.local.n[2],stride=na*nb;
  for(int mode=0;mode<=np/2;mode++)for(int v=0;v<4;v++){
   long double initial=0,remaining=0,boundary=0,focus_boundary=0;
   double maximum=0;int maximum_i=-1,maximum_j=-1;
   for(int slot=0;slot<np;slot++)if((slot<=np/2?slot:slot-np/2)==mode)
    for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
     long double bmode=0,rmode=0;
     for(int k=0;k<np;k++){
      int q=4*(i+na*j+stride*k)+v;double w=s.derivatives.forward[slot*np+k];
      bmode+=(long double)w*rhs[q];rmode+=(long double)w*(rhs[q]-action[q]);
     }
     initial+=bmode*bmode;remaining+=rmode*rmode;
     if(std::abs(rmode)>maximum){maximum=std::abs(rmode);maximum_i=i;maximum_j=j;}
     if(HISPID_FOCUS_TAU&&mode<=4)for(int side=0;side<hispid::focus_rank;side++)
      if(i+na*j==hispid::focus_row(side,na,nb))focus_boundary+=rmode*rmode;
     if(HISPID_AXIS_TAU&&mode>=5&&(i==0||j==0||j==nb-1))boundary+=rmode*rmode;
    }
   std::fprintf(stderr,"HiSpID residual mode=%d component=%d initial_l2=%.17g remaining_l2=%.17g tau_l2=%.17g focus_l2=%.17g maximum=%.17g maximum_i=%d maximum_j=%d\n",
      mode,v,double(std::sqrt(initial)),double(std::sqrt(remaining)),double(std::sqrt(boundary)),double(std::sqrt(focus_boundary)),maximum,maximum_i,maximum_j);
  }
  std::fprintf(stderr,"HiSpID action probe iterations=%d recurrence=%.17g true=%.17g relative=%.17g\n",
    result.iterations,result.recurrence_residual,result.true_residual,result.relative_residual);
  std::fflush(stderr);
 }
 const char*probe_setting=std::getenv("HISPID_PROBE_MODAL_FACTORS");
 if(probe_setting&&std::strcmp(probe_setting,"1")==0){
  std::fprintf(stderr,"HiSpID linear probe iterations=%d status=%d recurrence=%.17g true=%.17g relative=%.17g\n",
    result.iterations,status,result.recurrence_residual,result.true_residual,result.relative_residual);
  std::fflush(stderr);
 }
 if(status!=PK_SUCCESS&&!(status==PK_CALLBACK&&!hispid::last_error.empty()))
   hispid::last_error=PK_status_string(status);
 return status==PK_SUCCESS;
}
// Private compatibility facade used by the existing manufactured controls.
bool gmres(HiSpID_Data&s,const Sparse&M,const std::vector<double>&rhs,std::vector<double>&x,double rtol,bool eager=false){
 return linear_solve(s,M,rhs,x,rtol,PK_GMRES,eager);
}
void update_diag(HiSpID_Data&s,const std::vector<double>&res){
 for(int k=0;k<4;k++)s.diag.scaled_linf[k]=s.diag.unscaled_linf[k]=0;
 for(int p=0;p<s.npt;p++)for(int k=0;k<4;k++){
  if(!std::isfinite(res[4*p+k])){s.diag.scaled_linf[k]=s.diag.unscaled_linf[k]=INFINITY;continue;}
  s.diag.scaled_linf[k]=std::max(s.diag.scaled_linf[k],std::abs(res[4*p+k]));
  s.diag.unscaled_linf[k]=std::max(s.diag.unscaled_linf[k],std::abs(res[4*p+k]/row_weight(s,p)));
 }
}
void to_local(const HiSpID_Data&s,const double*x,double*y){
 for(int i=0;i<3;i++){y[i]=0;for(int j=0;j<3;j++)y[i]+=s.frame[j][i]*(x[j]-s.origin[j]);}
}
void make_coefficients(HiSpID_Data&s){
 if(s.coefficients_valid)return;
 auto start=std::chrono::steady_clock::now();s.coefficients.resize(s.ntotal);
#ifdef PUNCTURES_KOKKOS
 if(s.device&&s.setup.geometry_execution){s.device->make_coefficients(s.values.data(),s.coefficients.data());}
 else
#endif
 {
 std::vector<double>temporary(s.ntotal);
 s.derivatives.raw.along(0,s.derivatives.coefficient[0],4,s.values.data(),temporary.data());
 s.derivatives.raw.along(1,s.derivatives.coefficient[1],4,temporary.data(),s.coefficients.data());
 }
 s.setup.coefficient_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 s.coefficients_valid=true;
}
void polynomial_basis(int N,double x,std::vector<double>&value,std::vector<double>&first){
 value.resize(N);first.resize(N);value[0]=1;first[0]=0;value[1]=x;first[1]=1;
 for(int i=2;i<N;i++){value[i]=2*x*value[i-1]-value[i-2];first[i]=2*value[i-1]+2*x*first[i-1]-first[i-2];}
 value[0]=.5;
}
// Return P, P_t, P_eta for one orthonormal Fourier mode.
void evaluate_mode(const HiSpID_Data&s,int k,double t,double eta,double P[4][3]){
 const int na=s.local.n[0],nb=s.local.n[1];std::vector<double>Ta,Da,Tb,Db;
 const double lambda=hispid::AxisDerivatives::radial_stretch,den=lambda+(1-lambda)*t,sigma=hispid::AxisDerivatives::radial_inverse(t),kappa=hispid::AxisDerivatives::angular_stretch,T=std::tanh(kappa),zeta=std::atanh(eta*T)/kappa;
 double ts,tss;hispid::AxisDerivatives::radial_forward(sigma,ts,tss);
 polynomial_basis(na,2*sigma-1,Ta,Da);polynomial_basis(nb,zeta,Tb,Db);std::memset(P,0,12*sizeof(double));
 for(int j=0;j<nb;j++)for(int i=0;i<na;i++)for(int v=0;v<4;v++){
  const double c=s.coefficients[v+4*(i+na*(j+nb*k))];P[v][0]+=c*Ta[i]*Tb[j];P[v][1]+=((hispid::AxisDerivatives::infinity_power==1&&!hispid::AxisDerivatives::radial_polynomial&&!hispid::AxisDerivatives::dual_endpoint_map)?2*lambda/(den*den):2/ts)*c*Da[i]*Tb[j];P[v][2]+=(T/(kappa*(1-T*T*eta*eta)))*c*Ta[i]*Db[j];
 }
}
using ModeValues=std::array<std::array<double,3>,4>;
struct MeridionalValues {
 double sx2,sr2,sx,cx,sr,cr,a,B,rho;
 std::vector<ModeValues> modes;
};
// A sphere centered at the map origin and aligned to its axis has identical
// meridional coordinates around each polar ring. Only the Fourier sum and
// seed geometry vary with phi; the expensive two-dimensional P sums do not.
MeridionalValues sphere_ring(HiSpID_Data&s,double axial,double rho){
 make_coefficients(s);MeridionalValues out;out.rho=rho;
 const double w0=axial/s.b,q=rho/s.b,t=.5*(w0*w0+q*q-1),root=std::hypot(t,q);
 out.sx2=t<0?q*q/(root-t):root+t;out.sr2=t>0?q*q/(root+t):root-t;
 out.sx=std::sqrt(out.sx2);out.cx=std::sqrt(1+out.sx2);out.sr=std::sqrt(out.sr2);
 out.cr=w0/out.cx;out.a=out.sx/(out.cx+1);out.B=-out.cr/(1+out.sr);
 const double eta=-2*out.B/(1+out.B*out.B);out.modes.resize(s.local.n[2]);
 for(int k=0;k<s.local.n[2];k++){
  double values[4][3];evaluate_mode(s,k,out.a*out.a,eta,values);
  for(int v=0;v<4;v++)for(int d=0;d<3;d++)out.modes[k][v][d]=values[v][d];
 }
 return out;
}
void evaluate_coefficients(const HiSpID_Data&s,double a,double B,double phi,double V[4][4],double stable_sinR=-1,const MeridionalValues*ring=nullptr){
 const int np=s.local.n[2],half=np/2;const double t=a*a,eta=-2*B/(1+B*B),sn=stable_sinR>=0?stable_sinR:(1-B*B)/(1+B*B);
 std::memset(V,0,16*sizeof(double));
 for(int k=0;k<np;k++){
  const int m=k<=half?k:k-half,r=hispid::AxisDerivatives::exponent(m);double P[4][3];
  if(ring){for(int v=0;v<4;v++)for(int d=0;d<3;d++)P[v][d]=ring->modes[k][v][d];}
  else evaluate_mode(s,k,t,eta,P);
  const double norm=std::sqrt((m==0||m==half?1.:2.)/np),F=norm*(k<=half?std::cos(m*phi):std::sin(m*phi)),DF=norm*m*(k<=half?-std::sin(m*phi):std::cos(m*phi));
  const double S=(1+a)*std::pow(a*sn,r),Sa=.5*(std::pow(a,r)+ (r? r*(1+a)*std::pow(a,r-1):0))*std::pow(sn,r);
  const double etaB=-2*sn/(1+B*B),snB=-4*B/std::pow(1+B*B,2);
  const double Sb=r?(1+a)*std::pow(a,r)*r*std::pow(sn,r-1)*snB:0;
  for(int v=0;v<4;v++){
   V[v][0]+=S*P[v][0]*F;V[v][1]+=(Sa*P[v][0]+S*a*P[v][1])*F;
   V[v][2]+=(Sb*P[v][0]+S*etaB*P[v][2])*F;V[v][3]+=S*P[v][0]*DF;
  }
 }
}
void sample_fields(HiSpID_Data&s,const double*x,Fields&f,const MeridionalValues*ring=nullptr){
 make_coefficients(s);
 const double rho=ring?ring->rho:std::hypot(x[1],x[2]);
 /* Exact Cartesian first-derivative limits of the C2 modal interpolant.
  * Only m0 and m1 contribute. A tiny axis band avoids losing transverse
  * coordinates to rounding during coordinate inversion. */
 if(rho<1e-10*s.b){
  double A=-1,B,axial,transverse;const bool inner=std::abs(x[0])<s.b;
  if(inner){
   const double R=std::acos(std::clamp(x[0]/s.b,-1.0,1.0));
   B=std::tan(R/2-Piq);
   axial=-(1+B*B)/(2*s.b*std::sin(R));transverse=1/(s.b*std::sin(R));
  }else{
   if(std::abs(x[0])==s.b)throw std::runtime_error("puncture map focus is excluded from correction sampling");
   const double X=std::acosh(std::abs(x[0])/s.b);
   A=2*std::tanh(X/2)-1;B=x[0]>0?-1:1;
   const double a=.5*(A+1),sign=x[0]>0?1:-1;
   axial=sign*(1-a*a)/(s.b*std::sinh(X));transverse=sign/(s.b*std::sinh(X));
  }
  double cosine[4][4],sine[4][4];evaluate_coefficients(s,.5*(A+1),B,0,cosine);evaluate_coefficients(s,.5*(A+1),B,Pih,sine);
  for(int v=0;v<4;v++){
   f[v][0]=(A-1)*cosine[v][0];
   f[v][1]=inner?(A-1)*cosine[v][2]*axial:(cosine[v][0]+(A-1)*cosine[v][1])*axial;
   f[v][2]=inner?(cosine[v][0]+(A-1)*cosine[v][1])*transverse:(A-1)*cosine[v][2]*transverse;
   f[v][3]=inner?(sine[v][0]+(A-1)*sine[v][1])*transverse:(A-1)*sine[v][2]*transverse;
   f[v][0]+=f[v][2]*x[1]+f[v][3]*x[2];
  }
  return;
 }
 /* Stable inverse of x=b cosh(X) cos(R), rho=b sinh(X) sin(R).
  * Rationalize the small root rather than subtracting nearly equal
  * puncture distances (which loses R close to the axis). */
 const double w0=x[0]/s.b,q=rho/s.b,t=.5*(w0*w0+q*q-1),root=std::hypot(t,q);
 const double sx2=ring?ring->sx2:(t<0?q*q/(root-t):root+t);
 const double sr2=ring?ring->sr2:(t>0?q*q/(root+t):root-t);
 const double sx=ring?ring->sx:std::sqrt(sx2),cx=ring?ring->cx:std::sqrt(1+sx2),sr=ring?ring->sr:std::sqrt(sr2),cr=ring?ring->cr:w0/cx,a=ring?ring->a:sx/(cx+1),B=ring?ring->B:-cr/(1+sr);
 double phi=std::atan2(x[2],x[1]);if(phi<0)phi+=2*Pi;
 double V[4][4];evaluate_coefficients(s,a,B,phi,V,sr,ring);
 const double denominator=s.b*(sx2+sr2),AX=1-a*a,BR=.5*(1+B*B),co=x[1]/rho,si=x[2]/rho;
 for(int v=0;v<4;v++){
  const double uA=V[v][0]-2*(1-a)*V[v][1],uB=-2*(1-a)*V[v][2],up=-2*(1-a)*V[v][3];
  f[v][0]=-2*(1-a)*V[v][0];f[v][1]=(uA*AX*sx*cr-uB*BR*cx*sr)/denominator;
  const double radial=(uA*AX*cx*sr+uB*BR*sx*cr)/denominator;
  f[v][2]=radial*co-up*si/rho;f[v][3]=radial*si+up*co/rho;
 }

}
}
extern "C" {
const char *HiSpID_preconditioner_description(const HiSpID_Data*s){
 if(!s)return nullptr;
#ifdef PUNCTURES_KOKKOS
 if(constraint_blocks_enabled(*s,s->resolved_options.krylov))return s->half_turn?"constraint_block_upper_triangular_half_turn":"constraint_block_upper_triangular";
#endif
 if(HISPID_CURVED_MODAL_PRECONDITIONER)return "modal_curved_mixed_equilibrated";if(HISPID_LEFT_GMRES&&s->seed_family==HISPID_SEED_TRUMPET_R0_M&&s->resolved_options.krylov==PK_GMRES)return "modal_left_original_norm_confirmed";return relaxation_sweeps(*s)>1?"modal_two_sweep_omega0.4":"modal";
}
const char *HiSpID_residual_scaling(){
#if HISPID_AXIS_TAU
 return HISPID_ROW_POWER==3?"sin3_alpha_beta_with_axis_tau":"sin6_alpha_beta_with_axis_tau";
#endif
 if constexpr(HISPID_ROW_POWER==3)return HISPID_INFINITY_EQUILIBRATION?"sin3_alpha_beta_times_one_minus_t_pow_minus6":"sin3_alpha_beta";
 return HISPID_INFINITY_EQUILIBRATION?"sin6_alpha_beta_times_one_minus_t_pow_minus6":"sin6_alpha_beta";
}
const char *HiSpID_unknown_parameterization(){
 if constexpr(HISPID_FOCUS_TAU){
  static const auto value=[](){std::array<char,192> value{};
   std::snprintf(value.data(),value.size(),hispid::focus_rank==2?"modal_P_C2tauC4_focus_map_v9_r%.17g_k%.17g":"modal_P_C2tauC4_focus2_map_v10_r%.17g_k%.17g",hispid::AxisDerivatives::radial_stretch,hispid::AxisDerivatives::angular_stretch);return value;}();return value.data();}

 if constexpr(hispid::AxisDerivatives::dual_endpoint_map){
  static const auto value=[](){std::array<char,192> value{};
   std::snprintf(value.data(),value.size(),"modal_P_C2tauC4_beta14_map_v8_r%.17g_k%.17g",hispid::AxisDerivatives::radial_stretch,hispid::AxisDerivatives::angular_stretch);return value;}();return value.data();}

 if constexpr(hispid::AxisDerivatives::radial_polynomial){
  static const auto value=[](){std::array<char,192> value{};
   std::snprintf(value.data(),value.size(),"modal_P_C2tauC4_poly2_map_v7_r%.17g_k%.17g",hispid::AxisDerivatives::radial_stretch,hispid::AxisDerivatives::angular_stretch);return value;}();
  return value.data();}

 if constexpr(hispid::AxisDerivatives::infinity_power==2){
  static const auto value=[](){std::array<char,192> value{};
   std::snprintf(value.data(),value.size(),"modal_P_C2tauC4_inf2_map_v6_r%.17g_k%.17g",hispid::AxisDerivatives::radial_stretch,hispid::AxisDerivatives::angular_stretch);return value;}();
  static_assert(HISPID_AXIS_TAU||hispid::AxisDerivatives::infinity_power==1,"quadratic infinity map requires axis tau");return value.data();}

#if HISPID_AXIS_TAU
 static const auto tau=[](){std::array<char,128> value{};
  std::snprintf(value.data(),value.size(),"modal_P_C2tauC4_map_v5_r%.17g_k%.17g",hispid::AxisDerivatives::radial_stretch,hispid::AxisDerivatives::angular_stretch);return value;}();
 return tau.data();
#endif
 // Experimental C4 data must never be interpreted as the default C2 basis.
 // Consumers compare the full token, so default C2 readers reject C4 data.
 if constexpr(hispid::AxisDerivatives::regularity_cap==6){
  static const auto c4=[](){std::array<char,128> value{};
   std::snprintf(value.data(),value.size(),"modal_P_C4prolate_map_v4_r%.17g_k%.17g",hispid::AxisDerivatives::radial_stretch,hispid::AxisDerivatives::angular_stretch);return value;}();
  return c4.data();
 }
 if constexpr(hispid::AxisDerivatives::radial_stretch==.2&&hispid::AxisDerivatives::angular_stretch==2.)return "modal_P_C2prolate_mapped_v2";
 static const auto identifier=[](){std::array<char,128> value{};
  std::snprintf(value.data(),value.size(),"modal_P_C2prolate_map_v3_r%.17g_k%.17g",hispid::AxisDerivatives::radial_stretch,hispid::AxisDerivatives::angular_stretch);return value;}();
 return identifier.data();
}
int HiSpID_collocation_maps(double*out){
 if(!out)return -1;out[0]=hispid::AxisDerivatives::radial_stretch;out[1]=hispid::AxisDerivatives::angular_stretch;return 0;
}
static HiSpID_Data *create_context(const HiSpID_Config*c,bool sampler_only,int execution=0,int geometry_execution=0,int family=HISPID_SEED_QI,bool operator_only=false,const HiSpID_FreeData*free_data=nullptr,const HiSpID_Wave22*wave=nullptr){
#ifdef PUNCTURES_KOKKOS
 puncture::ExecutionLock execution_guard(puncture::execution_mutex(),std::defer_lock);if(execution)execution_guard.lock();
#endif
 hispid::last_error.clear();
 if(family!=HISPID_SEED_QI && family!=HISPID_SEED_TRUMPET_R0_M && family!=HISPID_SEED_LIU){hispid::last_error="invalid seed family";return nullptr;}
 if(sampler_only && execution){hispid::last_error="sampling-only contexts use CPU execution";return nullptr;}
 if(execution<0||execution>1){hispid::last_error="invalid execution backend";return nullptr;}
 if(geometry_execution<0||geometry_execution>1||(geometry_execution&&(!execution||sampler_only))){hispid::last_error="execution geometry requires a Kokkos solving context";return nullptr;}
#ifndef PUNCTURES_KOKKOS
 if(execution){hispid::last_error="library was built without Kokkos";return nullptr;}
#endif
 if(!c||!hispid::valid(*c,sampler_only,operator_only)){hispid::last_error="invalid HiSpID configuration";return nullptr;}
 if(free_data||wave){
  if(geometry_execution||family==HISPID_SEED_TRUMPET_R0_M||c->far_radius!=0||
     (c->hole[0].mass>0)==(c->hole[1].mass>0)){
   hispid::last_error="free data require one stationary QI/Liu hole, far_radius=0, host geometry";return nullptr;
  }
  for(const auto&h:c->hole)for(double v:h.velocity)if(v!=0){hispid::last_error="free data require stationary Kerr gauge";return nullptr;}
  if(free_data&&(free_data->struct_size!=sizeof(*free_data)||!free_data->evaluate)){hispid::last_error="invalid free-data provider";return nullptr;}
  if(wave){
   const int h=c->hole[0].mass>0?0:1;const auto&hole=c->hole[h];
   const double a=hispid::point_norm<double>(hole.spin)/hole.mass;
   const double throat=family==HISPID_SEED_LIU?(hole.mass+std::sqrt(hole.mass*hole.mass-a*a))/4:std::sqrt(hole.mass*hole.mass-a*a)/2;
   if(!hispid::valid_wave22(*wave)||wave->cutoff_min<=std::max(throat,c->inner_max[h])){
    hispid::last_error="invalid wave22 parameters or cutoff intersects Kerr throat/interior transition";return nullptr;
   }
  }
 }
 if(HISPID_FOCUS_TAU&&!sampler_only&&!operator_only){
  // A nonzero far-reference correction W can itself be non-smooth at a
  // boosted/trumpet puncture. These homogeneous rows constrain P alone.
  if(c->far_radius!=0){hispid::last_error="focus regularity requires far_radius=0";return nullptr;}
  if(!c->inner_flatten){hispid::last_error="focus regularity requires a flat interior plateau";return nullptr;}
  for(int h=0;h<2;h++)if(c->hole[h].mass>0&&!(c->inner_min[h]>0)){
   hispid::last_error="focus regularity requires positive inner_min for each active hole";return nullptr;
  }
 }
 const bool galerkin=HISPID_OVERINTEGRATED_GALERKIN&&!operator_only&&!sampler_only&&family==HISPID_SEED_TRUMPET_R0_M;
 if(galerkin&&geometry_execution){hispid::last_error="Galerkin currently requires the qualified host seed geometry";return nullptr;}
 HiSpID_Config qc=*c;qc.n[0]=(3*c->n[0]+1)/2;qc.n[1]=(3*c->n[1]+1)/2;
 long double allocation_bound=operator_only?3192.L*c->n[0]*c->n[1]*c->n[2]:hispid::allocation_bound(*c,sampler_only);
 const long double curved_bytes=HISPID_CURVED_MODAL_PRECONDITIONER&&!sampler_only&&!operator_only?(6.L*sizeof(hispid::ModalSymbol)+128.L)*(c->n[2]/2+1)*c->n[0]*c->n[1]:0;
 allocation_bound+=curved_bytes;
 if(galerkin&&!execution){
  allocation_bound+=3264.L*qc.n[0]*qc.n[1]*qc.n[2]+16.L*(c->n[0]*qc.n[0]+c->n[1]*qc.n[1]);
  if(allocation_bound>(long double)c->memory_limit_mib*1024*1024){hispid::last_error="Galerkin aggregate allocation bound exceeds memory_limit_mib";return nullptr;}
 }
#ifdef PUNCTURES_KOKKOS
 if(execution&&!sampler_only){
  try{puncture::initialize();}catch(const std::exception&e){hispid::last_error=e.what();return nullptr;}
  const long double points=(long double)c->n[0]*c->n[1]*c->n[2],radial=(long double)c->n[0]*c->n[0]*c->n[1]*2*(c->n[2]/2+1);
  const long double cpu=points*(3192+16*HISPID_STABLE_SCALAR_SOURCE+4*(5*12+48)+32*30+16.L/c->n[2])+radial*16+2.L*(c->n[2]/2+1)*c->n[1]*24*c->n[0];
  long double device=points*(sizeof(Cached)+8.L*(40+40+8+64+15+8+4*(2*c->krylov_restart+5)))+radial*24;
  if(HISPID_CONSTRAINT_BLOCK_PRECONDITIONER&&!operator_only&&!galerkin){
   const long double tables=12.L*c->n[0]*c->n[1]*2*(c->n[2]/2+1)*(HISPID_CURVED_MODAL_PRECONDITIONER?3:1);
   const long double fixed=points*(sizeof(Cached)+8.L*(40+40+8+64+15+12))+tables;
   // LU/inverse/transfer coexist only during prepare(); LU is freed before
   // either Krylov basis is allocated. The solve holds two matrix banks,
   // the outer basis, and at most one packed inner-block basis at a time.
   // Reserve the additional LGMRES directions for any supported outer method.
   const long double krylov=points*(32.L*(2*c->krylov_restart+16)+24.L*(ConstraintBlockInverse::restart+6+2*ConstraintBlockInverse::augmentation+1)+ConstraintBlockInverse::workspace_bytes_per_point);
   device=fixed+radial*16+std::max(radial*8,krylov);
   if constexpr(Kokkos::SpaceAccessibility<puncture::Exec,Kokkos::HostSpace>::accessible)
    device+=points*24.L*(ConstraintBlockInverse::restart+ConstraintBlockInverse::augmentation);
  }
  // Retained CPU data, device data (or extra host copies for OpenMP), lazy
  // maximum Krylov basis, packed factor mirrors and LU/inverse overlap.
  long double small=0;for(int axis=0;axis<3;axis++)small+=16.L*c->n[axis]*c->n[axis];small+=32.L*c->n[2]*c->n[2]+160.L*c->n[0]*c->n[1];
  allocation_bound=cpu+device+radial*16+small+curved_bytes+
   (HISPID_CONSTRAINT_BLOCK_PRECONDITIONER&&!operator_only&&!galerkin&&!Kokkos::SpaceAccessibility<puncture::Exec,Kokkos::HostSpace>::accessible?points*24.L*(ConstraintBlockInverse::restart+ConstraintBlockInverse::augmentation):0)+2*(hispid::polar_border_bytes(c->n[0],c->n[1],c->n[2])+hispid::focus_border_bytes(c->n[0],c->n[1],c->n[2]))+HISPID_AXIS_TAU*64.L*(2*c->n[0]+c->n[1])*c->n[2];
  // Kokkos owns its Krylov basis on the execution space; there is no
  // simultaneous host PK_solve basis. Keep all host outer vectors above.
  if(HISPID_CURVED_MODAL_PRECONDITIONER&&!operator_only)device+=(2.L*sizeof(hispid::ModalSymbol)+32.L)*(c->n[2]/2+1)*c->n[0]*c->n[1];
  const long double raw_device=points*(sizeof(Cached)+8.L*(40+40+8+64+15+8));
  if(operator_only){allocation_bound=points*3192+raw_device+small;device=raw_device;}
  if(!operator_only&&family==HISPID_SEED_TRUMPET_R0_M&&(HISPID_DAMPED_MODAL_RELAXATION||HISPID_LEFT_GMRES))device+=(HISPID_LEFT_GMRES?128.L:96.L)*points;
  if(galerkin){
   const long double qp=(long double)qc.n[0]*qc.n[1]*qc.n[2];
   const long double matrices=32.L*(c->n[0]*qc.n[0]+c->n[1]*qc.n[1]);
   const long double qdevice=qp*(sizeof(Cached)+8.L*(40+40+8+64+15+8));
   // Quadrature contexts never factor a preconditioner or allocate a
   // Krylov basis. On CUDA the host geometry cache is released after upload.
   long double retained=qp*32+qdevice+qp*32+matrices;
   if constexpr(Kokkos::SpaceAccessibility<puncture::Exec,Kokkos::HostSpace>::accessible)retained+=qp*3192;
   const long double setup_peak=points*(3192+sizeof(Cached)+8.L*(40+40+8+64+15+8))+qp*3192+qdevice+small+matrices;
   allocation_bound=std::max(allocation_bound+retained,setup_peak);
   device+=qdevice+qp*32+matrices;
  }
  if(allocation_bound>(long double)c->memory_limit_mib*1024*1024){hispid::last_error="Kokkos aggregate allocation bound exceeds memory_limit_mib";return nullptr;}
#ifdef KOKKOS_ENABLE_CUDA
  if constexpr(std::is_same_v<puncture::Exec,Kokkos::Cuda>){
   auto free=Puncture_execution_device_free_bytes();
   if(!free||device+small+(hispid::polar_border_bytes(c->n[0],c->n[1],c->n[2])+hispid::focus_border_bytes(c->n[0],c->n[1],c->n[2]))+HISPID_AXIS_TAU*32.L*(2*c->n[0]+c->n[1])*c->n[2]+512.L*1024*1024>free){hispid::last_error="Kokkos device allocation bound exceeds available GPU memory or query failed";return nullptr;}
  }
#endif
 }
#endif
 if(!sampler_only&&!operator_only&&family==HISPID_SEED_TRUMPET_R0_M&&(HISPID_DAMPED_MODAL_RELAXATION||HISPID_LEFT_GMRES))allocation_bound+=192.L*c->n[0]*c->n[1]*c->n[2];
 if(!sampler_only)allocation_bound+=HISPID_SPECTRAL_RADIAL_PRECONDITIONER*(16.L*c->n[0]*c->n[0]+24.L*c->n[0]);
 if(allocation_bound>(long double)c->memory_limit_mib*1024*1024){hispid::last_error="radial preconditioner allocation bound exceeds memory_limit_mib";return nullptr;}
 HiSpID_Data*s=nullptr;try{
  s=new HiSpID_Data;s->allocation_bound=allocation_bound;s->seed_family=family;s->config=s->local=*c;s->sampler_only=sampler_only;s->operator_only=operator_only;s->setup.geometry_execution=geometry_execution;
  if(free_data)s->free_data=*free_data;if(wave){s->wave=*wave;s->has_wave=true;}
  if(HISPID_HALF_TURN_SOLVE&&!sampler_only){
   if(family!=HISPID_SEED_LIU||s->free_data.evaluate||((c->hole[0].mass>0)==(c->hole[1].mass>0))||
      c->hole[0].center[0]!=c->hole[1].center[0]||c->hole[0].center[1]!=c->hole[1].center[1])
    throw std::runtime_error("half-turn solve requires a single Liu hole and a lab-z aligned grid without callbacks");
   for(const auto&h:c->hole)if(h.mass>0){
    if(h.spin[0]!=0||h.spin[1]!=0||h.velocity[0]!=0||h.velocity[1]!=0||h.velocity[2]!=0)
     throw std::runtime_error("half-turn solve requires an unboosted lab-z aligned spin");
   }
   s->half_turn=true;
  }
  if(geometry_execution)s->setup.scalar_digits=std::numeric_limits<double>::digits;
  double sep[3],len=0;for(int i=0;i<3;i++){s->origin[i]=.5*(c->hole[0].center[i]+c->hole[1].center[i]);sep[i]=c->hole[0].center[i]-c->hole[1].center[i];len+=sep[i]*sep[i];}
  len=std::sqrt(len);
  if(!std::isfinite(len))throw std::runtime_error("nonfinite puncture separation");
  if(len<1e-12)throw std::runtime_error("puncture map centers must be separated, including an inactive focus");
  s->b=len/2;for(int i=0;i<3;i++)s->frame[i][0]=sep[i]/len;
  int least=0;for(int i=1;i<3;i++)if(std::abs(s->frame[i][0])<std::abs(s->frame[least][0]))least=i;
  double norm=std::sqrt(1-s->frame[least][0]*s->frame[least][0]);
  for(int i=0;i<3;i++)s->frame[i][1]=((i==least?1.0:0.0)-s->frame[i][0]*s->frame[least][0])/norm;
  for(int i=0;i<3;i++)s->frame[i][2]=s->frame[(i+1)%3][0]*s->frame[(i+2)%3][1]-s->frame[(i+2)%3][0]*s->frame[(i+1)%3][1];
  for(int h=0;h<2;h++){
   to_local(*s,c->hole[h].center,s->local.hole[h].center);
   for(int i=0;i<3;i++){s->local.hole[h].spin[i]=s->local.hole[h].velocity[i]=0;for(int j=0;j<3;j++){
    s->local.hole[h].spin[i]+=s->frame[j][i]*c->hole[h].spin[j];s->local.hole[h].velocity[i]+=s->frame[j][i]*c->hole[h].velocity[j];}}
  }
  s->npt=c->n[0]*c->n[1]*c->n[2];s->ntotal=4*s->npt;s->values.resize(s->ntotal);
  s->diag.npoints=s->npt;
  if(sampler_only){s->derivatives.initialize(c->n,true);return s;}
#ifdef PUNCTURES_KOKKOS
  if(geometry_execution){
   s->device=std::make_unique<HiKokkos>(c->n,s->derivatives,s->geometry,s->b,&s->local,family,operator_only);
   s->device->row_weights(s->row_weights);
   s->setup.spectral_seconds=s->device->spectral_setup_seconds;s->setup.geometry_seconds=s->device->geometry_setup_seconds;
   return s;
  }
#endif
  s->geometry.resize(s->npt);if(!execution)allocate_derivs(&s->work,s->ntotal);
  auto setup_start=std::chrono::steady_clock::now();
  s->derivatives.initialize(c->n);
  s->setup.spectral_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-setup_start).count();
  auto build_geometry=[&](int p){
   const int i=p%c->n[0],j=(p/c->n[0])%c->n[1],k=p/(c->n[0]*c->n[1]);
   Fields f{};double x[3];transform(*s,i,j,k,f,x);hispid::Background bg;context_background(*s,x,bg);
   cache(bg,s->geometry[p]);
   double sn=std::sin(Pih*(2*i+1)/c->n[0])*std::sin(Pih*(2*j+1)/c->n[1]);s->geometry[p].weight=std::pow(sn,HISPID_ROW_POWER);
   if constexpr(HISPID_INFINITY_EQUILIBRATION){
    const double a=.5*(s->derivatives.coordinate[0][i]+1);
    s->geometry[p].weight/=std::pow(1-a*a,6);
   }
  };
  setup_start=std::chrono::steady_clock::now();
#ifdef PUNCTURES_KOKKOS
  if(execution){
   puncture::initialize();std::mutex lock;std::string error;
   Kokkos::parallel_for("host seed geometry",Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0,s->npt),[&](int p){try{build_geometry(p);}catch(const std::exception&e){std::lock_guard<std::mutex>guard(lock);if(error.empty())error=e.what();}});
   Kokkos::DefaultHostExecutionSpace().fence();if(!error.empty())throw std::runtime_error(error);
   s->setup.geometry_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-setup_start).count();
   s->device=std::make_unique<HiKokkos>(c->n,s->derivatives,s->geometry,s->b,nullptr,family,operator_only);
   s->setup.spectral_seconds+=s->device->spectral_setup_seconds;s->setup.geometry_seconds+=s->device->geometry_setup_seconds;
   s->row_weights.resize(c->n[0]*c->n[1]);for(size_t p=0;p<s->row_weights.size();p++)s->row_weights[p]=s->geometry[p].weight;
   if constexpr(!Kokkos::SpaceAccessibility<puncture::Exec,Kokkos::HostSpace>::accessible){s->geometry.clear();s->geometry.shrink_to_fit();}
  }else
#endif
  {for(int p=0;p<s->npt;p++)build_geometry(p);
   s->setup.geometry_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-setup_start).count();}
  s->diag.npoints=s->npt;
  if(galerkin){
   s->quadrature.reset(create_context(&qc,false,execution,geometry_execution,family,true));
   if(!s->quadrature)throw std::runtime_error("Galerkin quadrature: "+hispid::last_error);
   s->transfer=std::make_unique<hispid::GalerkinTransfer>(c->n,qc.n);
#ifdef PUNCTURES_KOKKOS
   if(execution)s->device_transfer=std::make_unique<DeviceGalerkinTransfer>(*s->transfer);
#endif
   s->setup.geometry_seconds+=s->quadrature->setup.geometry_seconds;s->setup.spectral_seconds+=s->quadrature->setup.spectral_seconds;
  }
 }catch(const std::exception&e){hispid::last_error=e.what();HiSpID_destroy(s);return nullptr;}return s;
}
HiSpID_Data *HiSpID_create_with_free_data(const HiSpID_Config*c,int family,int execution,int sampler_only,const HiSpID_FreeData*f){
 if(!f||(sampler_only!=0&&sampler_only!=1)){hispid::last_error="invalid free-data constructor arguments";return nullptr;}
 return create_context(c,sampler_only,execution,0,family,false,f);
}
HiSpID_Data *HiSpID_create_with_wave22(const HiSpID_Config*c,int family,int execution,int sampler_only,const HiSpID_Wave22*w){
 if(!w||(sampler_only!=0&&sampler_only!=1)){hispid::last_error="invalid wave22 constructor arguments";return nullptr;}
 return create_context(c,sampler_only,execution,0,family,false,nullptr,w);
}
int HiSpID_wave22_sample(const HiSpID_Wave22*w,const double*center,int count,const double*xyz,double*h,double*ht){
 if(!w||!center||count<0||!xyz||!h||!ht||!hispid::valid_wave22(*w))return -1;
 for(int i=0;i<3;i++)if(!std::isfinite(center[i]))return -1;
 for(int p=0;p<count;p++){
  for(int i=0;i<3;i++)if(!std::isfinite(xyz[3*p+i]))return -1;
  hispid::Jet a[3][3],b[3][3];hispid::wave22(*w,center,xyz+3*p,a,b);
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){h[9*p+3*i+j]=a[i][j].v;ht[9*p+3*i+j]=b[i][j].v;}
 }return 0;
}
HiSpID_Data *HiSpID_create_with_seed_family(const HiSpID_Config*c,int family,int execution,int geometry_execution,int sampler_only){
 if(sampler_only!=0 && sampler_only!=1){hispid::last_error="invalid sampler mode";return nullptr;}
 return create_context(c,sampler_only!=0,execution,geometry_execution,family);
}
int HiSpID_seed_family(const HiSpID_Data*s){return s?s->seed_family:-1;}
HiSpID_Data *HiSpID_create(const HiSpID_Config*c){return create_context(c,false);}
HiSpID_Data *HiSpID_create_with_execution(const HiSpID_Config*c,int execution){return create_context(c,false,execution);}
HiSpID_Data *HiSpID_create_with_geometry(const HiSpID_Config*c,int execution,int geometry_execution){return create_context(c,false,execution,geometry_execution);}
int HiSpID_setup_statistics(const HiSpID_Data*s,HiSpID_SetupStatistics*out){
 if(!s||!out||out->struct_size!=sizeof(HiSpID_SetupStatistics))return -1;*out=s->setup;return 0;
}
HiSpID_Data *HiSpID_create_sampler(const HiSpID_Config*c){return create_context(c,true);}
static bool sampling_context(HiSpID_Data*s){
 if(s&&s->sampler_only){hispid::last_error="sampling-only context cannot evaluate or solve collocation equations";return true;}
 return false;
}
static int solve_context(HiSpID_Data*s,double fixed_forcing,int method=PK_GMRES){
#ifdef PUNCTURES_KOKKOS
 auto execution_guard=execution_lock(s);
#endif

 if(sampling_context(s))return -1;
 if(!s)return -1;
 if(s->operator_only){hispid::last_error="operator-only context cannot solve";return -1;}
 hispid::last_error.clear();auto start=std::chrono::steady_clock::now();
 if(method==PK_LGMRES){
  // Three appended Arnoldi/search pairs, three saved corrections and scratch.
  // One MiB additionally bounds scalar banks for the allowed restart <=200.
  const long double extra=80.L*s->ntotal+1024.L*1024;
  if(s->allocation_bound+extra>(long double)s->local.memory_limit_mib*1024*1024){
   hispid::last_error="LGMRES aggregate allocation bound exceeds memory_limit_mib";return -1;
  }
#ifdef KOKKOS_ENABLE_CUDA
  if(s->device){
   if constexpr(std::is_same_v<puncture::Exec,Kokkos::Cuda>){
    auto free=Puncture_execution_device_free_bytes();
    if(!free||extra+512.L*1024*1024>free){hispid::last_error="LGMRES workspace exceeds available GPU memory";return -1;}
   }
  }
#endif
 }
 s->resolved_options={int(sizeof(HiSpID_SolveOptions)),method,fixed_forcing};
 s->diag={};s->diag.npoints=s->npt;s->coefficients_valid=false;
 s->jvp_applications=s->preconditioner_applications=0;s->linear_history.clear();
 const char*trace_setting=std::getenv("HISPID_TRACE_NEWTON");
 const bool trace=trace_setting&&std::strcmp(trace_setting,"1")==0;
 std::vector<double>r(s->ntotal),rhs(s->ntotal),step,trial(s->ntotal),rt(s->ntotal);
 // The normalized vector FD matrices depend only on this solve's fixed
 // map/grid/mode. Move their factors between Newton steps; scalar factors
 // are rebuilt from the current nonlinear potential at every step.
 std::vector<ModalBlock>vector_cache;
 try{
  if(s->half_turn)std::fprintf(stderr,"HiSpID half-turn unknown sector enabled; original full RHS/JVP/residual gates retained\n");
  project_solve_sector(*s,s->values.data(),true);
  residual(*s,s->values.data(),r.data());
  for(int it=0;it<=s->local.max_newton;it++){
   if(s->half_turn){long double total=0,defect=0;const int N=s->local.n[0]*s->local.n[1],np=s->local.n[2];
    for(int k=0;k<np/2;k++)for(int q=0;q<4*N;q++){double a=r[q+4*N*k],b=r[q+4*N*(k+np/2)],e=a-(q%4>=2?-b:b);defect+=(long double)e*e;total+=(long double)a*a+(long double)b*b;}
    std::fprintf(stderr,"HiSpID full residual half-turn it=%d forbidden_l2=%.17g original_l2=%.17g\n",it,std::sqrt(double(defect/2)),std::sqrt(double(total)));
   }
   update_diag(*s,r);double err=norminf(r);if(err<=s->local.tolerance){s->diag.converged=1;break;}
   if(it==s->local.max_newton)break;
   s->diag.newton_iterations++;
   Sparse M=preconditioner(*s,&vector_cache);for(int i=0;i<s->ntotal;i++)rhs[i]=-r[i];
   project_solve_sector(*s,rhs.data(),false);
   double forcing=fixed_forcing>0?fixed_forcing:std::min(.05,std::max(1e-5,std::sqrt(err)));
   const int krylov_before=s->diag.krylov_iterations;
   const bool linear_ok=linear_solve(*s,M,rhs,step,forcing,method);M.retain_vectors(vector_cache);
   s->linear_history.push_back({double(it),forcing,s->last_gmres_relative,double(s->diag.krylov_iterations-krylov_before)});
   if(!linear_ok)break;
   bool accepted=false;double old=norm2v(r);
   std::vector<double>trace_action;
   double trace_action_norm=0;
   if(trace){
    trace_action.resize(s->ntotal);jvp(*s,step.data(),trace_action.data());
    trace_action_norm=norm2v(trace_action);
    for(int i=0;i<s->ntotal;i++)rt[i]=r[i]+trace_action[i];
    std::fprintf(stderr,"HiSpID Newton trace it=%d residual_l2=%.17g step_linf=%.17g linear_relative=%.17g slope=%.17g\n",
      it,old,norminf(step),norm2v(rt)/old,dot(r,trace_action)/(old*old));
    std::fflush(stderr);
   }
   for(int backtrack=0;backtrack<=HISPID_NEWTON_BACKTRACKS;backtrack++){
    const double damping=std::ldexp(1.0,-backtrack);
    for(int i=0;i<s->ntotal;i++)trial[i]=s->values[i]+damping*step[i];
    residual(*s,trial.data(),rt.data());
    if(trace){
     auto&observed=s->quadrature?*s->quadrature:*s;
     double remainder=0,min_psi=std::numeric_limits<double>::infinity();int nonfinite=0;
     for(int i=0;i<s->ntotal;i++){
      const double q=rt[i]-r[i]-damping*trace_action[i];remainder+=q*q;
      if(!std::isfinite(rt[i]))nonfinite++;
     }
#ifdef PUNCTURES_KOKKOS
     if(s->device){
      auto c=observed.device->geometry;auto b=observed.device->base;
      Kokkos::parallel_reduce("Newton trace minimum psi",puncture::Range(0,observed.npt),
        KOKKOS_LAMBDA(int p,double&v){const double psi=c(p).psi+b(40*p);if(psi<v)v=psi;},Kokkos::Min<double>(min_psi));
     }else
#endif
     for(int p=0;p<observed.npt;p++)min_psi=std::min(min_psi,observed.geometry[p].psi+observed.basefields[p][0][0]);
     std::fprintf(stderr,"HiSpID Newton trial it=%d damping=%.17g residual_ratio=%.17g directional_remainder=%.17g min_psi=%.17g nonfinite=%d\n",
       it,damping,norm2v(rt)/old,std::sqrt(remainder)/(damping*trace_action_norm),min_psi,nonfinite);
     std::fflush(stderr);
    }
    if(norm2v(rt)<old*(1-1e-4*damping)&&std::isfinite(norminf(rt))){s->values=trial;r=rt;accepted=true;break;}
   }
   if(!accepted){hispid::last_error="Newton line search failed";residual(*s,s->values.data(),r.data());break;}
  }
 }catch(const std::exception&e){hispid::last_error=e.what();}
 try{residual(*s,s->values.data(),r.data());update_diag(*s,r);}
 catch(const std::exception&e){hispid::last_error=e.what();s->diag.converged=0;}
 s->diag.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 if(!s->diag.converged&&hispid::last_error.empty())hispid::last_error="Newton iteration limit";
 return s->diag.converged?0:1;
}
int HiSpID_solve(HiSpID_Data*s){return solve_context(s,0);}
int HiSpID_set_axisymmetric(HiSpID_Data*s,int enabled){
 if(!s||s->sampler_only||(enabled!=0&&enabled!=1))return -1;
 if(enabled&&s->half_turn){hispid::last_error="half-turn solve retains full equations and cannot enable axisymmetric equation projection";return -1;}
 if(enabled&&(s->has_wave||s->free_data.evaluate)){hispid::last_error="custom free data require the full 3D solve";return -1;}
 if(enabled)for(const auto&h:s->local.hole){
  if(h.mass<=0)continue;
  for(double spin:h.spin)if(spin!=0){hispid::last_error="axisymmetric no-swirl solve requires zero spins";return -1;}
  if(h.velocity[1]!=0||h.velocity[2]!=0||h.center[1]!=0||h.center[2]!=0){
   hispid::last_error="axisymmetric solve requires exactly coaxial local centers and boosts";return -1;
  }
 }
 s->axisymmetric=enabled;return 0;
}
int HiSpID_solve_with_forcing(HiSpID_Data*s,double rtol){
 if(!std::isfinite(rtol)||rtol<=0||rtol>=1){hispid::last_error="invalid fixed relative linear tolerance";return -1;}
 return solve_context(s,rtol);
}
void HiSpID_default_solve_options(HiSpID_SolveOptions*out){if(out)*out={int(sizeof(*out)),PK_GMRES,0};}
int HiSpID_solve_with_options(HiSpID_Data*s,const HiSpID_SolveOptions*o){
 if(!o||o->struct_size!=sizeof(*o)||(o->krylov!=PK_GMRES&&o->krylov!=PK_BICGSTAB&&o->krylov!=PK_LGMRES)||
    !std::isfinite(o->linear_rtol)||o->linear_rtol<0||o->linear_rtol>=1){hispid::last_error="invalid linear solve options";return -1;}
 return solve_context(s,o->linear_rtol,o->krylov);
}
int HiSpID_work_statistics(const HiSpID_Data*s,int*out){
 if(!s||!out)return -1;out[0]=s->jvp_applications;out[1]=s->preconditioner_applications;return 0;
}
int HiSpID_resolved_solve_options(const HiSpID_Data*s,HiSpID_SolveOptions*out){
 if(!s||!out||out->struct_size!=sizeof(*out)||s->resolved_options.struct_size!=sizeof(*out))return -1;
 *out=s->resolved_options;return 0;
}
int HiSpID_linear_history(const HiSpID_Data*s,int capacity,double*out){
 if(!s||capacity<0)return -1;int n=s->linear_history.size();
 if(!out)return capacity==0?n:-1;
 if(capacity==0||capacity<n)return -1;
 if(out)for(int i=0;i<n;i++)for(int j=0;j<4;j++)out[4*i+j]=s->linear_history[i][j];
 return n;
}
int HiSpID_diagnostics(const HiSpID_Data*s,HiSpID_Diagnostics*d){if(!s||!d)return -1;*d=s->diag;return 0;}
int HiSpID_get_unknowns(const HiSpID_Data*s,double*v,int n){if(!s||!v||n!=s->ntotal)return -1;std::copy(s->values.begin(),s->values.end(),v);return 0;}
int HiSpID_set_unknowns(HiSpID_Data*s,const double*v,int n){if(!s||!v||n!=s->ntotal)return -1;for(int i=0;i<n;i++)if(!std::isfinite(v[i]))return -1;std::copy(v,v+n,s->values.begin());s->coefficients_valid=false;s->diag.converged=0;return 0;}
int HiSpID_residual(HiSpID_Data*s,const double*v,double*r){
#ifdef PUNCTURES_KOKKOS
 auto execution_guard=execution_lock(s);
#endif

 if(sampling_context(s))return -1;
 if(!s||!v||!r)return -1;
 try{residual(*s,v,r);}catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
int HiSpID_jvp(HiSpID_Data*s,const double*v,const double*d,double*r){
#ifdef PUNCTURES_KOKKOS
 auto execution_guard=execution_lock(s);
#endif

 if(sampling_context(s))return -1;
 if(!s||!v||!d||!r)return -1;
 if(s->quadrature){
  try{std::vector<double>scratch(s->ntotal);residual(*s,v,scratch.data());jvp(*s,d,r);}catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
 }
#ifdef PUNCTURES_KOKKOS
 if(s->device&&s->setup.geometry_execution){
  try{s->device->fields(v,true);jvp(*s,d,r);}catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
 }
#endif
 try{fields(*s,v,s->basefields);jvp(*s,d,r);}catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
int HiSpID_equation_samples(HiSpID_Data*s,double*xyz,double*g,double*psi,double*HM){
#ifdef PUNCTURES_KOKKOS
 auto execution_guard=execution_lock(s);
#endif

 if(sampling_context(s))return -1;
 if(!s||!xyz||!g||!psi||!HM)return -1;
 try{
#ifdef PUNCTURES_KOKKOS
  if(s->device&&s->setup.geometry_execution){
   s->device->fields(s->values.data(),true);auto c=s->device->geometry;
   auto base=s->device->base,position=s->device->op.positions,trig=s->device->op.trig;
   const int na=s->local.n[0],nb=s->local.n[1];Kokkos::Array<double,3>origin;Kokkos::Array<double,9>frame;
   for(int a=0;a<3;a++){origin[a]=s->origin[a];for(int b=0;b<3;b++)frame[3*a+b]=s->frame[a][b];}
   puncture::View samples("collocation diagnostic export",9ull*s->npt);
   Kokkos::parallel_for("collocation diagnostic export",puncture::Range(0,s->npt),KOKKOS_LAMBDA(int p){
    int row=p%(na*nb),k=p/(na*nb);double local[3]={position(2*row),position(2*row+1)*trig(2*k),position(2*row+1)*trig(2*k+1)};
    Fields u{};for(int v=0;v<4;v++)for(int d=0;d<10;d++)u[v][d]=base(40*p+10*v+d);
    double factor=c(p).psi+u[0][0],res[4];eval(c(p),u,res);res[0]*=-8/std::pow(factor,5);for(int d=1;d<4;d++)res[d]/=std::pow(factor,10);
    for(int a=0;a<3;a++){double x=origin[a],momentum=0;for(int b=0;b<3;b++){x+=frame[3*a+b]*local[b];momentum+=frame[3*a+b]*res[b+1];}samples(9*p+a)=x;samples(9*p+6+a)=momentum;}
    samples(9*p+3)=c(p).g;samples(9*p+4)=factor;samples(9*p+5)=res[0];
   });
   std::vector<double>host(9ull*s->npt);puncture::download(samples,host.data());
   for(int p=0;p<s->npt;p++){for(int a=0;a<3;a++)xyz[3*p+a]=host[9*p+a];g[p]=host[9*p+3];psi[p]=host[9*p+4];for(int a=0;a<4;a++)HM[4*p+a]=host[9*p+5+a];}
   return 0;
  }
#endif
  fields(*s,s->values.data(),s->basefields);
  std::vector<Cached>mirror;const Cached*geometry=s->geometry.data();
#ifdef PUNCTURES_KOKKOS
  if(s->device&&s->geometry.empty()){mirror.resize(s->npt);using H=Kokkos::View<Cached*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;Kokkos::deep_copy(H(mirror.data(),mirror.size()),s->device->geometry);geometry=mirror.data();}
#endif
  for(int k=0;k<s->local.n[2];k++)for(int j=0;j<s->local.n[1];j++)for(int i=0;i<s->local.n[0];i++){
   const int p=pindex(*s,i,j,k);const Cached&c=geometry[p];
   psi[p]=c.psi+s->basefields[p][0][0];g[p]=c.g;
   eval(c,s->basefields[p],HM+4*p);
   HM[4*p]*=-8/std::pow(psi[p],5);
   for(int d=1;d<4;d++)HM[4*p+d]/=std::pow(psi[p],10);
   double momentum[3]={};for(int a=0;a<3;a++)for(int b=0;b<3;b++)momentum[a]+=s->frame[a][b]*HM[4*p+b+1];
   for(int a=0;a<3;a++)HM[4*p+a+1]=momentum[a];
   Fields zero{};double x[3];transform(*s,i,j,k,zero,x);
   for(int a=0;a<3;a++){
    xyz[3*p+a]=s->origin[a];for(int b=0;b<3;b++)xyz[3*p+a]+=s->frame[a][b]*x[b];
   }
  }
 }catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
static int sample_with_derivatives(HiSpID_Data*s,int count,const double*xyz,HiSpID_Point*out,double*dgamma,const MeridionalValues*ring=nullptr,const double*local_xyz=nullptr){
 if(!s||!xyz||!out||count<0)return -1;
 try{for(int p=0;p<count;p++){
  double x[3];if(local_xyz)std::copy(local_xyz+3*p,local_xyz+3*p+3,x);else to_local(*s,xyz+3*p,x);
  hispid::Background bg;context_background(*s,x,bg);Fields f{};sample_fields(*s,x,f,ring);
  f[0][0]+=bg.far_correction.v;
  double psi=bg.psi.v+f[0][0];if(!(psi>0))throw std::runtime_error("nonpositive solved conformal factor");
  double L[3][3],gam[3][3],K[3][3],A[3][3];L_and_div(bg.metric,bg.inv,bg.C,f,L,nullptr);
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   A[i][j]=bg.M[i][j].v+L[i][j];gam[i][j]=std::pow(psi,4)*(double)bg.metric[i][j].v;
   K[i][j]=A[i][j]/(psi*psi)+gam[i][j]*(double)bg.K.v/3;
  }
  std::memset(out+p,0,sizeof(*out));out[p].psi=psi;out[p].mean_curvature=bg.K.v;out[p].attenuation=bg.g.v;
  out[p].correction[0]=f[0][0];for(int i=0;i<3;i++)for(int a=0;a<3;a++)out[p].correction[i+1]+=s->frame[i][a]*f[a+1][0];
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int a=0;a<3;a++)for(int b=0;b<3;b++){
   double q=s->frame[i][a]*s->frame[j][b];out[p].gamma[3*i+j]+=q*gam[a][b];out[p].Kij[3*i+j]+=q*K[a][b];
   out[p].conformal_metric[3*i+j]+=q*(double)bg.metric[a][b].v;out[p].Atilde[3*i+j]+=q*A[a][b];
  }
  if(dgamma){
   std::fill(dgamma+27*p,dgamma+27*(p+1),0.0);
   for(int d=0;d<3;d++)for(int i=0;i<3;i++)for(int j=0;j<3;j++)
    for(int c=0;c<3;c++)for(int a=0;a<3;a++)for(int b=0;b<3;b++){
     const double dpsi=bg.psi.d[c+1]+f[0][c+1]+bg.far_correction.d[c+1];
     const double grad=4*std::pow(psi,3)*dpsi*(double)bg.metric[a][b].v
                        +std::pow(psi,4)*(double)bg.metric[a][b].d[c+1];
     dgamma[27*p+9*d+3*i+j]+=s->frame[d][c]*s->frame[i][a]*s->frame[j][b]*grad;
    }
  }
 }}catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
int HiSpID_sample(HiSpID_Data*s,int count,const double*xyz,HiSpID_Point*out){
 return sample_with_derivatives(s,count,xyz,out,nullptr);
}
int HiSpID_sample_with_derivatives(HiSpID_Data*s,int count,const double*xyz,HiSpID_Point*out,double*dgamma){
 if(!dgamma)return -1;
 return sample_with_derivatives(s,count,xyz,out,dgamma);
}
int HiSpID_charges(HiSpID_Data*s,const double*center,double radius,int nt,int np,double*out){
 if(!s||!center||!out||radius<=0||nt<4||np<8)return -1;
 try{
 std::fill(out,out+7,0.0);
 const bool centered=center[0]==s->origin[0]&&center[1]==s->origin[1]&&center[2]==s->origin[2];
 /* Gauss-Legendre cos(theta) and uniform phi. */
 for(int a=0;a<nt;a++){
  double mu=std::cos(Pi*(a+.75)/(nt+.5)),pp=0;
  for(int it=0;it<20;it++){
   double p0=1,p1=mu;for(int l=2;l<=nt;l++){double p=((2*l-1)*mu*p1-(l-1)*p0)/l;p0=p1;p1=p;}
   pp=nt*(mu*p1-p0)/(mu*mu-1);double step=p1/pp;mu-=step;if(std::abs(step)<1e-15)break;
  }double weight=2/((1-mu*mu)*pp*pp)*(2*Pi/np)*radius*radius;
  MeridionalValues ring;if(centered)ring=sphere_ring(*s,radius*mu,radius*std::sqrt(1-mu*mu));
  for(int b=0;b<np;b++){
   // Align the integration polar axis with the prolate separation axis.
   // This keeps high meridional polynomial degrees out of the azimuthal
   // quadrature. Normals/tensors and returned E,P,J remain in the lab frame.
   double ph=2*Pi*(b+.5)/np,local_normal[3]={mu,std::sqrt(1-mu*mu)*std::cos(ph),std::sqrt(1-mu*mu)*std::sin(ph)},n[3]={};
   for(int i=0;i<3;i++)for(int j=0;j<3;j++)n[i]+=s->frame[i][j]*local_normal[j];
   double x[3];for(int i=0;i<3;i++)x[i]=center[i]+radius*n[i];HiSpID_Point p;double dg[3][9];
   // Reuse the tested analytic physical metric gradient. This evaluates
   // the retained modal polynomial once rather than at13 FD stencil points.
   // Independent physical-FD charge comparisons remain validation controls.
   double local_x[3];for(int i=0;i<3;i++)local_x[i]=radius*local_normal[i];
   if(sample_with_derivatives(s,1,x,&p,&dg[0][0],centered?&ring:nullptr,centered?local_x:nullptr))return -2;
   double E=0,P[3]={};for(int i=0;i<3;i++)for(int j=0;j<3;j++){
    E+=n[i]*(dg[j][3*i+j]-dg[i][3*j+j]);
    P[i]+=(p.Kij[3*i+j]-p.mean_curvature*p.gamma[3*i+j])*n[j];
   }
   out[0]+=weight*E/(16*Pi);for(int i=0;i<3;i++){
    out[i+1]+=weight*P[i]/(8*Pi);out[i+4]+=weight*radius*(n[(i+1)%3]*P[(i+2)%3]-n[(i+2)%3]*P[(i+1)%3])/(8*Pi);
   }
  }
 }return 0;
 }catch(const std::exception&e){hispid::last_error=e.what();return -2;}
}
void HiSpID_destroy(HiSpID_Data*s){
#ifdef PUNCTURES_KOKKOS
 auto execution_guard=execution_lock(s);
#endif
if(!s)return;if(s->work)free_derivs(s->work);delete s;}
}
