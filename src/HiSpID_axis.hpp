#ifndef HISPID_AXIS_HPP
#define HISPID_AXIS_HPP
#include "HiSpID_spectral.hpp"
#include "HiSpID_jets.hpp"
#include <algorithm>
#include <type_traits>
#ifndef HISPID_RADIAL_STRETCH
#define HISPID_RADIAL_STRETCH .2
#endif
#ifndef HISPID_ANGULAR_STRETCH
#define HISPID_ANGULAR_STRETCH 2.
#endif

#ifndef HISPID_DUAL_ENDPOINT_MAP
#define HISPID_DUAL_ENDPOINT_MAP 0
#endif
#ifndef HISPID_RADIAL_POLYNOMIAL
#define HISPID_RADIAL_POLYNOMIAL 0
#endif
#ifndef HISPID_INFINITY_POWER
#define HISPID_INFINITY_POWER 1
#endif
#ifndef HISPID_REGULARITY_CAP
#define HISPID_REGULARITY_CAP 4
#endif
namespace hispid {
// Regular prolate modes: V_m=(1+a)[a sin(R)]^r P_m(a²,cos(R)).
// Default r=m for m<=4, then 3/4 for odd/even m. This preserves every smooth
// Cartesian mode and guarantees C2 without dividing by arbitrarily high
// powers near the axes. Unknowns are the nodal values of modal P; Fourier
// modes use the orthonormal real basis below. No near-axis Fourier division.
struct AxisDerivatives {
 std::array<int,3> n{};
 SpectralDerivatives raw;
 std::array<std::vector<double>,2> coordinate,first_map,second_map;
 std::array<std::vector<double>,2> coefficient;
 std::vector<double> forward,inverse,inverse_phi,inverse_phi2;
 mutable std::array<std::vector<double>,6> scratch;
 static constexpr double radial_stretch=HISPID_RADIAL_STRETCH,angular_stretch=HISPID_ANGULAR_STRETCH;
 // The rational map is monotone for every positive stretch. Values above
 // one retain an exterior radial scale with a nearby inactive focus.
 static_assert(radial_stretch>=.001&&radial_stretch<=64,"radial map stretch must be in[.001,64]");
 static_assert(angular_stretch>=.1&&angular_stretch<=6,"angular map stretch must be in[.1,6]");
 static constexpr bool dual_endpoint_map=HISPID_DUAL_ENDPOINT_MAP;
 static_assert(!dual_endpoint_map||radial_stretch<=.5,"dual endpoint slope must be in[.001,.5]");
 static constexpr bool radial_polynomial=HISPID_RADIAL_POLYNOMIAL;
 static_assert(!radial_polynomial||radial_stretch<=1,"polynomial radial map stretch must be <=1");
 static constexpr int infinity_power=HISPID_INFINITY_POWER;
 static_assert(infinity_power==1||infinity_power==2,"infinity power must be1or2");
 // Positive polynomial density. F is the Beta(15,3) CDF; its
 // binomial-tail form avoids cancellation at both endpoints. t_s(0)=l,
 // t_s(1)=2. No extra vanishing power at compactified infinity is introduced.
 template<class R> HISPID_GEOMETRY_INLINE static R dual_value(R s,R l){
  const R b=(R(2)-l)/R(5),c=R(1)-l-b,q=R(1)-s;
  return l*s+b*std::pow(s,5)+c*std::pow(s,15)*(R(136)*q*q+R(17)*s*q+s*s);
 }
 template<class R> HISPID_GEOMETRY_INLINE static void dual_slopes(R s,R l,R&ts,R&tss){
  const R b=(R(2)-l)/R(5),c=R(1)-l-b,q=R(1)-s;
  ts=l+R(5)*b*std::pow(s,4)+c*R(2040)*std::pow(s,14)*q*q;
  tss=R(20)*b*std::pow(s,3)+c*R(2040)*(R(14)*std::pow(s,13)*q*q-R(2)*std::pow(s,14)*q);
 }
 template<class R> HISPID_GEOMETRY_INLINE static R dual_inverse(R t){
  if(t<=R(0))return R(0);if(t>=R(1))return R(1);
  R lo=0,hi=1,s=std::sqrt(t),l=R(radial_stretch);
  for(int k=0;k<64;k++){
   R value=dual_value(s,l),ts,tss;if(value==t)return s;
   dual_slopes(s,l,ts,tss);if(value<t)lo=s;else hi=s;
   R next=s-(value-t)/ts;if(next==s)return s;
   if(!(next>lo&&next<hi))next=(lo+hi)/R(2);
   if(next==s)return s;s=next;
  }return s;
 }
 // t=a^2; (1-t) ~ (1-s)^p. Squaring the infinity distance raises
 // differentiability of polyhomogeneous terms without changing free data.
 HISPID_GEOMETRY_INLINE static double radial_t(double s){
  if constexpr(dual_endpoint_map)return dual_value(s,radial_stretch);
  if constexpr(radial_polynomial)return radial_stretch*s+(1-radial_stretch)*s*s;
  if constexpr(infinity_power==1)return radial_stretch*s/(1-(1-radial_stretch)*s);
  return radial_stretch*s/(radial_stretch*s+(1-s)*(1-s));
 }
 template<class T> HISPID_GEOMETRY_INLINE static T radial_inverse(T t){
  if constexpr(dual_endpoint_map){
   if constexpr(std::is_arithmetic_v<T>)return dual_inverse(t);
   else{auto s=dual_inverse(t.v),ts=s,tss=s;dual_slopes(s,decltype(s)(radial_stretch),ts,tss);return unary(t,s,1/ts,-tss/(ts*ts*ts));}
  }
  if constexpr(radial_polynomial){using std::sqrt;return T(2)*t/(T(radial_stretch)+sqrt(T(radial_stretch*radial_stretch)+T(4*(1-radial_stretch))*t));}
  if constexpr(infinity_power==1)return t/(T(radial_stretch)+T(1-radial_stretch)*t);
  T q=T(radial_stretch)*(T(1)-t);using std::sqrt;
  return T(2)*t/(T(2)*t+q+sqrt(q*(T(4)*t+q)));
 }
 HISPID_GEOMETRY_INLINE static void radial_forward(double s,double&ts,double&tss){
  const double l=radial_stretch;
  if constexpr(dual_endpoint_map)dual_slopes(s,l,ts,tss);
  else if constexpr(radial_polynomial){ts=l+2*(1-l)*s;tss=2*(1-l);}
  else if constexpr(infinity_power==1){double d=1-(1-l)*s;ts=l/(d*d);tss=2*(1-l)*l/(d*d*d);}
  else{double d=l*s+(1-s)*(1-s),ds=l-2+2*s;ts=l*(1-s*s)/(d*d);tss=-2*l*s/(d*d)-2*l*(1-s*s)*ds/(d*d*d);}
 }
 HISPID_GEOMETRY_INLINE static double angular_eta(double z){return std::tanh(angular_stretch*z)/std::tanh(angular_stretch);}
 static constexpr int regularity_cap=HISPID_REGULARITY_CAP;
 static_assert(regularity_cap==4||regularity_cap==6,"regularity cap must be4or6");
 HISPID_GEOMETRY_INLINE static int exponent(int m){return m<=regularity_cap?m:(m%2?regularity_cap-1:regularity_cap);}
 HISPID_GEOMETRY_INLINE static double radial(double a,int r){return (1+a)*std::pow(a,r);}
 HISPID_GEOMETRY_INLINE static double angular(double eta,int r){double q=1-eta*eta;return std::pow(q>0?q:0.,.5*r);}
 HISPID_GEOMETRY_INLINE static void node(int axis,int N,int i,double out[3]){
  const double pi=std::acos(-1.0),z=-std::cos(pi*(i+.5)/N);
  if(axis==0){const double t=radial_t(.5*(1+z)),a=std::sqrt(t),den=radial_stretch+(1-radial_stretch)*t;
   out[0]=2*a-1;
   if constexpr(infinity_power==1&&!radial_polynomial&&!dual_endpoint_map){out[1]=2*a*radial_stretch/(den*den);
    out[2]=radial_stretch/(den*den)-4*(1-radial_stretch)*t*radial_stretch/(den*den*den);
   }else{double ts,tss;radial_forward(.5*(1+z),ts,tss);out[1]=2*a/ts;out[2]=1/ts-2*t*tss/(ts*ts*ts);}
  }else{const double eta=angular_eta(z),B=-eta/(1+std::sqrt(1-eta*eta)),den=1+B*B,eb=-2*(1-B*B)/(den*den),ebb=4*B*(3-B*B)/(den*den*den),T=std::tanh(angular_stretch),Q=1-T*T*eta*eta;
   out[0]=B;out[1]=T*eb/(angular_stretch*Q);
   out[2]=T*ebb/(angular_stretch*Q)+2*T*T*T*eta*eb*eb/(angular_stretch*Q*Q);
  }
 }
 HISPID_GEOMETRY_INLINE static double coefficient_at(int N,int i,int j){return (2./N)*(i%2?-1:1)*std::cos(std::acos(-1.0)*i*(j+.5)/N);}
 HISPID_GEOMETRY_INLINE static void fourier(int N,int k,int mode,double out[3]){
  const int half=N/2,m=mode<=half?mode:mode-half;const bool cosine=mode<=half;
  const double phi=2*std::acos(-1.0)*k/N,normal=std::sqrt((m==0||m==half?1.:2.)/N);
  out[0]=normal*(cosine?std::cos(m*phi):std::sin(m*phi));
  out[1]=(m==0||m==half)?0:normal*m*(cosine?-std::sin(m*phi):std::cos(m*phi));
  out[2]=-m*m*out[0];
 }
 void initialize(const int*shape,bool sampling_only=false){
  std::copy(shape,shape+3,n.begin());
  for(int axis=0;axis<2;axis++){
   const int N=n[axis];coefficient[axis].resize(N*N);
   for(int i=0;i<N;i++)for(int j=0;j<N;j++)coefficient[axis][i*N+j]=coefficient_at(N,i,j);
   coordinate[axis].resize(N);first_map[axis].resize(N);second_map[axis].resize(N);
   for(int i=0;i<N;i++){
    double values[3];node(axis,N,i,values);coordinate[axis][i]=values[0];first_map[axis][i]=values[1];second_map[axis][i]=values[2];
   }
  }
  const int N=n[2];
  forward.resize(N*N);inverse.resize(N*N);inverse_phi.resize(N*N);inverse_phi2.resize(N*N);
  for(int k=0;k<N;k++)for(int mode=0;mode<N;mode++){
   double values[3];fourier(N,k,mode,values);forward[mode*N+k]=inverse[k*N+mode]=values[0];inverse_phi[k*N+mode]=values[1];inverse_phi2[k*N+mode]=values[2];
  }
  raw.n=n;
  if(sampling_only)return;
  raw.initialize(shape);
 }
 void apply(int nv,derivs*w)const{
  const int size=nv*n[0]*n[1]*n[2];for(auto&v:scratch)v.resize(size);
  std::copy(w->d0,w->d0+size,scratch[0].begin());
  raw.differentiate(0,raw.D[0],nv,w->d0,scratch[1].data());raw.differentiate(1,raw.D[1],nv,w->d0,scratch[2].data());
  raw.differentiate(0,raw.D2[0],nv,w->d0,scratch[3].data());raw.differentiate(1,raw.D2[1],nv,w->d0,scratch[4].data());
  raw.differentiate(0,raw.D[0],nv,scratch[2].data(),scratch[5].data());
  const int half=n[2]/2;
  for(int p=0;p<size/nv;p++){
   const int i=p%n[0],j=(p/n[0])%n[1],mode=p/(n[0]*n[1]),m=mode<=half?mode:mode-half,r=exponent(m);
   const double a=.5*(coordinate[0][i]+1),B=coordinate[1][j],eta=-2*B/(1+B*B),den=1-eta*eta;
   const double ma=first_map[0][i],mma=second_map[0][i],mb=first_map[1][j],mmb=second_map[1][j],S=radial(a,r)*angular(eta,r);
   const double la=.5*(1/(1+a)+r/a),lla=-.25*(1/((1+a)*(1+a))+r/(a*a));
   const double etaB=-2*(1-B*B)/std::pow(1+B*B,2),etaBB=4*B*(3-B*B)/std::pow(1+B*B,3);
   const double le=-r*eta/den,lle=-r*(1+eta*eta)/(den*den),lb=le*etaB,llb=lle*etaB*etaB+le*etaBB;
   for(int v=0;v<nv;v++){
    const int q=nv*p+v;const double P=scratch[0][q],Pa=scratch[1][q],Pb=scratch[2][q];
    w->d0[q]=S*P;w->d1[q]=S*(ma*Pa+la*P);w->d2[q]=S*(mb*Pb+lb*P);
    w->d11[q]=S*(ma*ma*scratch[3][q]+(mma+2*la*ma)*Pa+(la*la+lla)*P);
    w->d22[q]=S*(mb*mb*scratch[4][q]+(mmb+2*lb*mb)*Pb+(lb*lb+llb)*P);
    w->d12[q]=S*(ma*mb*scratch[5][q]+la*mb*Pb+lb*ma*Pa+la*lb*P);
   }
  }
  raw.along(2,inverse_phi,nv,w->d0,w->d3);raw.along(2,inverse_phi2,nv,w->d0,w->d33);
  raw.along(2,inverse_phi,nv,w->d1,w->d13);raw.along(2,inverse_phi,nv,w->d2,w->d23);
  for(double*value:{w->d0,w->d1,w->d2,w->d11,w->d22,w->d12}){
   raw.along(2,inverse,nv,value,scratch[0].data());std::copy(scratch[0].begin(),scratch[0].end(),value);
  }
 }
};
}
#endif
