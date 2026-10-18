// Frozen C1-only implementation, symbols renamed only for compatibility control.
#ifndef HISPID_FOCUS_REFERENCE_HPP
#define HISPID_FOCUS_REFERENCE_HPP
#include "HiSpID_axis.hpp"
#ifndef HISPID_FOCUS_TAU
#define HISPID_FOCUS_TAU 0
#endif
namespace hispid {
// Bounded Cartesian corrections in the flat interior are smooth. For m<=4,
// u_m/rho^m ~ (1-t)^(m+1) P_m, and smoothness at (t,eta)=(0,+/-1)
// requires (d_t-2*side*d_eta-(m+1))P_m=0. These replace just two rows.
inline std::vector<double> reference_focus_weights(int na,int nb,int mode){
 const int N=na*nb;std::vector<double>out(2*N);
 std::vector<long double>wa(na),da(na),wb(nb),db(nb);
 auto endpoint=[](int n,int side,std::vector<long double>&w,std::vector<long double>&d){
  const long double pi=std::acos(-1.L);std::vector<long double>dx(n);long double sum=0,s=0;
  for(int i=0;i<n;i++){long double th=pi*(i+.5L)/n;dx[i]=side<0?2*std::sin(th/2)*std::sin(th/2):-2*std::cos(th/2)*std::cos(th/2);w[i]=(i%2?-1.L:1.L)*std::sin(th)/(-dx[i]);sum+=w[i];}
  for(int i=0;i<n;i++){w[i]/=sum;s+=w[i]/dx[i];}
  for(int i=0;i<n;i++)d[i]=w[i]*(1/dx[i]-s);
 };
 endpoint(na,-1,wa,da);for(auto&x:da)x*=2/AxisDerivatives::radial_stretch;
 const long double k=AxisDerivatives::angular_stretch,T=std::tanh(k),eta_scale=T/(k*(1-T*T));
 for(int side=0;side<2;side++){
  int sign=side?1:-1;endpoint(nb,sign,wb,db);long double norm=0;
  for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
   long double q=da[i]*wb[j]-2*sign*wa[i]*eta_scale*db[j]-(mode+1)*wa[i]*wb[j];
   out[side*N+i+na*j]=q;norm+=std::abs(q);
  }
  for(int i=0;i<N;i++)out[side*N+i]/=norm;
 }return out;
}
inline std::vector<double> reference_focus_bank(int na,int nb,int np){
 std::vector<double>w;for(int m=0;m<=std::min(4,np/2);m++){auto row=reference_focus_weights(na,nb,m);w.insert(w.end(),row.begin(),row.end());}return w;
}
static inline long double reference_focus_border_bytes(int na,int nb,int np){
 return HISPID_FOCUS_TAU?128.L*na*nb*2*std::min(5,np/2+1)+256.L*np:0;
}
inline void reference_focus_differences(int na,int nb,int np,const double*w,const double*inv,const double*u,const double*r,std::vector<double>&delta){
 const int N=na*nb;delta.assign(8*np,0);
 for(int k=0;k<np;k++){int m=k<=np/2?k:k-np/2;if(m>4)continue;
  for(int side=0;side<2;side++)for(int v=0;v<4;v++){
   long double q=0;for(int p=0;p<N;p++)q+=(long double)w[(2*m+side)*N+p]*u[4*(p+N*k)+v];
   int row=side?na*(nb-1):0;for(int l=0;l<np;l++)q-=(long double)inv[l*np+k]*r[4*(row+N*l)+v];
   delta[4*(side*np+k)+v]=q;
  }
 }
}
inline void reference_focus_replace(int na,int nb,int np,const double*inv,const std::vector<double>&delta,double*r){
 const int N=na*nb;for(int side=0;side<2;side++)for(int k=0;k<np;k++)for(int v=0;v<4;v++){
  double q=0;for(int m=0;m<np;m++)q+=inv[k*np+m]*delta[4*(side*np+m)+v];
  r[4*((side?na*(nb-1):0)+N*k)+v]+=q;
 }
}
}
#endif
