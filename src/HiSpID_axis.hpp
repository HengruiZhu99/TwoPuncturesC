#ifndef HISPID_AXIS_HPP
#define HISPID_AXIS_HPP
#include "HiSpID_spectral.hpp"
#include <algorithm>
#ifndef HISPID_RADIAL_STRETCH
#define HISPID_RADIAL_STRETCH .2
#endif
#ifndef HISPID_ANGULAR_STRETCH
#define HISPID_ANGULAR_STRETCH 2.
#endif

namespace hispid {
// Regular prolate modes: V_m=(1+a)[a sin(R)]^r P_m(a²,cos(R)).
// r=m for m<=4, then 3/4 for odd/even m. This preserves every smooth
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
 static_assert(radial_stretch>=.001&&radial_stretch<=1,"radial map stretch must be in[.001,1]");
 static_assert(angular_stretch>=.1&&angular_stretch<=6,"angular map stretch must be in[.1,6]");
 static double radial_t(double s){return radial_stretch*s/(1-(1-radial_stretch)*s);}
 static double angular_eta(double z){return std::tanh(angular_stretch*z)/std::tanh(angular_stretch);}
 static int exponent(int m){return m<=4?m:(m%2?3:4);}
 static double radial(double a,int r){return (1+a)*std::pow(a,r);}
 static double angular(double eta,int r){return std::pow(std::max(0.,1-eta*eta),.5*r);}
 void initialize(const int*shape,bool sampling_only=false){
  std::copy(shape,shape+3,n.begin());const double pi=std::acos(-1.0);
  for(int axis=0;axis<2;axis++){
   const int N=n[axis];coefficient[axis].resize(N*N);
   for(int i=0;i<N;i++)for(int j=0;j<N;j++)coefficient[axis][i*N+j]=(2./N)*(i%2?-1:1)*std::cos(pi*i*(j+.5)/N);
   coordinate[axis].resize(N);first_map[axis].resize(N);second_map[axis].resize(N);
   for(int i=0;i<N;i++){
    const double z=-std::cos(pi*(i+.5)/N);
    if(axis==0){const double t=radial_t(.5*(1+z)),a=std::sqrt(t),den=radial_stretch+(1-radial_stretch)*t;
     coordinate[0][i]=2*a-1;first_map[0][i]=2*a*radial_stretch/(den*den);
     second_map[0][i]=radial_stretch/(den*den)-4*(1-radial_stretch)*t*radial_stretch/(den*den*den);}
    else{const double eta=angular_eta(z),B=-eta/(1+std::sqrt(1-eta*eta)),den=1+B*B,eb=-2*(1-B*B)/(den*den),ebb=4*B*(3-B*B)/(den*den*den),T=std::tanh(angular_stretch),Q=1-T*T*eta*eta;
     coordinate[1][i]=B;first_map[1][i]=T*eb/(angular_stretch*Q);
     second_map[1][i]=T*ebb/(angular_stretch*Q)+2*T*T*T*eta*eb*eb/(angular_stretch*Q*Q);}
   }
  }
  const int N=n[2],half=N/2;
  forward.resize(N*N);inverse.resize(N*N);inverse_phi.resize(N*N);inverse_phi2.resize(N*N);
  for(int k=0;k<N;k++)for(int mode=0;mode<N;mode++){
   const int m=mode<=half?mode:mode-half;const bool cosine=mode<=half;
   const double phi=2*pi*k/N,normal=std::sqrt((m==0||m==half?1.:2.)/N);
   const double value=normal*(cosine?std::cos(m*phi):std::sin(m*phi));
   forward[mode*N+k]=inverse[k*N+mode]=value;
   inverse_phi[k*N+mode]=(m==0||m==half)?0:normal*m*(cosine?-std::sin(m*phi):std::cos(m*phi));
   inverse_phi2[k*N+mode]=-m*m*value;
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
