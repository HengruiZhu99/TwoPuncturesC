#ifndef HISPID_FOCUS_HPP
#define HISPID_FOCUS_HPP
#include "HiSpID_axis.hpp"
#ifndef HISPID_FOCUS_TAU
#define HISPID_FOCUS_TAU 0
#endif
namespace hispid {
#ifndef HISPID_FOCUS_ORDER
#define HISPID_FOCUS_ORDER 1
#endif
static_assert(HISPID_FOCUS_ORDER==1||HISPID_FOCUS_ORDER==2,"focus order must be1or2");
constexpr int focus_rank=2*HISPID_FOCUS_ORDER;
HISPID_GEOMETRY_INLINE int focus_row(int q,int na,int nb){return q/2+(q%2?na*(nb-1):0);}
// Smooth Cartesian corrections in the flat interior obey two necessary
// compatibility identities for Q=(1-t)^(m+1)P at (t,eta)=(0,+/-1):
// Q_t=2*side*Q_eta, Q_tt=4*Q_etaeta+4*side*Q_eta.
// Select one or two collocation rows per focus, with L1-normalized identities.
inline std::vector<double> focus_weights(int na,int nb,int mode){
 const int N=na*nb;std::vector<double>out(focus_rank*N);
 std::vector<long double>wa(na),da(na),daa(na),wb(nb),db(nb),dbb(nb);
 auto endpoint=[](int n,int side,std::vector<long double>&w,std::vector<long double>&d,std::vector<long double>&dd){
  const long double pi=std::acos(-1.L);std::vector<long double>r(n);long double sum=0,s=0,s2=0;
  for(int i=0;i<n;i++){long double th=pi*(i+.5L)/n,dx=side<0?2*std::sin(th/2)*std::sin(th/2):-2*std::cos(th/2)*std::cos(th/2);r[i]=-1/dx;w[i]=(i%2?-1.L:1.L)*std::sin(th)*r[i];sum+=w[i];}
  for(int i=0;i<n;i++){w[i]/=sum;s+=w[i]*r[i];s2+=w[i]*r[i]*r[i];}
  for(int i=0;i<n;i++){d[i]=w[i]*(s-r[i]);dd[i]=2*w[i]*(s*s-s*r[i]+r[i]*r[i]-s2);}
 };
 endpoint(na,-1,wa,da,daa);const long double lam=AxisDerivatives::radial_stretch,zt=2/lam,ztt=-4*(1-lam)/(lam*lam);
 for(int i=0;i<na;i++){daa[i]=zt*zt*daa[i]+ztt*da[i];da[i]*=zt;}
 const long double k=AxisDerivatives::angular_stretch,T=std::tanh(k),eta_scale=T/(k*(1-T*T));
 for(int side=0;side<2;side++){
  int sign=side?1:-1;endpoint(nb,sign,wb,db,dbb);
  for(int j=0;j<nb;j++){dbb[j]=eta_scale*eta_scale*dbb[j]+2*sign*T*T*T/(k*(1-T*T)*(1-T*T))*db[j];db[j]*=eta_scale;}
  for(int order=0;order<focus_rank/2;order++){
   long double norm=0;const int row=2*order+side;
   for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
    long double q=order?daa[i]*wb[j]-4*wa[i]*dbb[j]-2*(mode+1)*da[i]*wb[j]-4*sign*wa[i]*db[j]+mode*(mode+1)*wa[i]*wb[j]:da[i]*wb[j]-2*sign*wa[i]*db[j]-(mode+1)*wa[i]*wb[j];
    out[row*N+i+na*j]=q;norm+=std::abs(q);
   }
   for(int i=0;i<N;i++)out[row*N+i]/=norm;
  }
 }return out;
}
inline std::vector<double> focus_bank(int na,int nb,int np){
 std::vector<double>w;for(int m=0;m<=std::min(4,np/2);m++){auto row=focus_weights(na,nb,m);w.insert(w.end(),row.begin(),row.end());}return w;
}
static inline long double focus_border_bytes(int na,int nb,int np){
 // Per low-mode block: weights+responses = 2*rank*N doubles. Up to
 // two component groups per mode; this bound covers both banks plus scratch.
 return HISPID_FOCUS_TAU?128.L*na*nb*2*std::min(5,np/2+1)+512.L*np:0;
}
inline void focus_differences(int na,int nb,int np,const double*w,const double*inv,const double*u,const double*r,std::vector<double>&delta){
 const int N=na*nb;delta.assign(4*focus_rank*np,0);
 for(int k=0;k<np;k++){int m=k<=np/2?k:k-np/2;if(m>4)continue;
  for(int side=0;side<focus_rank;side++)for(int v=0;v<4;v++){
   long double q=0;for(int p=0;p<N;p++)q+=(long double)w[(focus_rank*m+side)*N+p]*u[4*(p+N*k)+v];
   int row=focus_row(side,na,nb);for(int l=0;l<np;l++)q-=(long double)inv[l*np+k]*r[4*(row+N*l)+v];
   delta[4*(side*np+k)+v]=q;
  }
 }
}
inline void focus_replace(int na,int nb,int np,const double*inv,const std::vector<double>&delta,double*r){
 const int N=na*nb;for(int side=0;side<focus_rank;side++)for(int k=0;k<np;k++)for(int v=0;v<4;v++){
  double q=0;for(int m=0;m<np;m++)q+=inv[k*np+m]*delta[4*(side*np+m)+v];
  r[4*(focus_row(side,na,nb)+N*k)+v]+=q;
 }
}
}
#endif
