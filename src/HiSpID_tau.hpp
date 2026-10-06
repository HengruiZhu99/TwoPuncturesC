#ifndef HISPID_TAU_HPP
#define HISPID_TAU_HPP
#include "HiSpID_axis.hpp"
#ifndef HISPID_AXIS_TAU
#define HISPID_AXIS_TAU 0
#endif
namespace hispid {
// C2 modal variables; replace independent near-axis rows with endpoint P=0.
// Polar rows include corners; radial rows exclude them to avoid redundancy.
HISPID_GEOMETRY_INLINE int tau_count(int na,int nb){return 2*na+nb-2;}
HISPID_GEOMETRY_INLINE void tau_node(int q,int na,int nb,int&i,int&j){
 if(q<na){i=q;j=0;}else if(q<2*na){i=q-na;j=nb-1;}else{i=0;j=q-2*na+1;}
}
inline std::vector<double> tau_weights(int na,int nb){
 std::vector<double>w(2*nb+na);const double pi=std::acos(-1.);
 for(int side=0;side<3;side++){
  int n=side==2?na:nb,offset=side==2?2*nb:side*nb;double endpoint=side==1?1:-1,sum=0;
  for(int k=0;k<n;k++){double theta=pi*(k+.5)/n;w[offset+k]=(k%2?-1.:1.)*std::sin(theta)/(endpoint+std::cos(theta));sum+=w[offset+k];}
  for(int k=0;k<n;k++)w[offset+k]/=sum;
 }return w;
}
HISPID_GEOMETRY_INLINE double tau_delta(int q,int mode,int v,int na,int nb,int np,
 const double*w,const double*inverse,const double*u,const double*r){
 const int m=mode<=np/2?mode:mode-np/2;if(m<5)return 0;
 int i,j;tau_node(q,na,nb,i,j);double old=0,value=0;
 for(int k=0;k<np;k++)old+=inverse[k*np+mode]*r[4*(i+na*(j+nb*k))+v];
 if(j==0||j==nb-1){const int offset=j==0?0:nb;
  for(int l=0;l<nb;l++)value+=w[offset+l]*u[4*(i+na*(l+nb*mode))+v];
 }else for(int l=0;l<na;l++)value+=w[2*nb+l]*u[4*(l+na*(j+nb*mode))+v];
 return value-old;
}
inline void tau_apply(int na,int nb,int np,const double*w,const double*inverse,
 const double*u,double*r,std::vector<double>&delta){
 const int count=tau_count(na,nb);delta.resize(4*count*np);
 for(int q=0;q<count;q++)for(int m=0;m<np;m++)for(int v=0;v<4;v++)
  delta[4*(q*np+m)+v]=tau_delta(q,m,v,na,nb,np,w,inverse,u,r);
 for(int q=0;q<count;q++){int i,j;tau_node(q,na,nb,i,j);
  for(int k=0;k<np;k++)for(int v=0;v<4;v++){
   double correction=0;for(int m=0;m<np;m++)correction+=inverse[k*np+m]*delta[4*(q*np+m)+v];
   r[4*(i+na*(j+nb*k))+v]+=correction;
  }
 }
}
}
#endif
