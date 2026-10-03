#ifndef HISPID_JETS_HPP
#define HISPID_JETS_HPP
#include <cmath>
#include <stdexcept>
#ifdef PUNCTURES_KOKKOS
#include <Kokkos_Core.hpp>
#define HISPID_GEOMETRY_INLINE KOKKOS_INLINE_FUNCTION
#else
#define HISPID_GEOMETRY_INLINE inline
#endif
namespace hispid {
/* Second-order forward jets. The host reference retains long double;
 * device builders instantiate double. Shared analytical identities remove
 * throat/puncture cancellation where possible.
 * No finite differences in seed geometry. */
template<class Real> struct JetT {
  Real v, d[4], h[4][4];
  HISPID_GEOMETRY_INLINE JetT(Real x=0):v(x),d{},h{} {}
  HISPID_GEOMETRY_INLINE static JetT<Real> variable(Real v,int k){JetT<Real> a(v);a.d[k]=1;return a;}
};
using Jet=JetT<long double>;
template<class T> struct ScalarArgument {using type=T;};
template<class Real> HISPID_GEOMETRY_INLINE Real scalar_power(Real x,Real p){
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
 return ::pow(double(x),double(p));
#else
 return std::pow(x,p);
#endif
}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> operator+(const JetT<Real>&a,const JetT<Real>&b){
  JetT<Real> c(a.v+b.v);for(int i=0;i<4;i++){c.d[i]=a.d[i]+b.d[i];
    for(int j=0;j<4;j++)c.h[i][j]=a.h[i][j]+b.h[i][j];}return c;
}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> operator-(const JetT<Real>&a,const JetT<Real>&b){
  JetT<Real> c(a.v-b.v);for(int i=0;i<4;i++){c.d[i]=a.d[i]-b.d[i];
    for(int j=0;j<4;j++)c.h[i][j]=a.h[i][j]-b.h[i][j];}return c;
}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> operator-(const JetT<Real>&a){return JetT<Real>(0)-a;}
template<class Real> HISPID_GEOMETRY_INLINE void compensated_add(JetT<Real>&sum,JetT<Real>&correction,const JetT<Real>&term){
 JetT<Real> adjusted=term-correction,next=sum+adjusted;correction=(next-sum)-adjusted;sum=next;
}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> operator*(const JetT<Real>&a,const JetT<Real>&b){
  JetT<Real> c(a.v*b.v);for(int i=0;i<4;i++){c.d[i]=a.d[i]*b.v+a.v*b.d[i];
    for(int j=0;j<4;j++)c.h[i][j]=a.h[i][j]*b.v+a.d[i]*b.d[j]
      +a.d[j]*b.d[i]+a.v*b.h[i][j];}return c;
}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> unary(const JetT<Real>&a,typename ScalarArgument<Real>::type v,typename ScalarArgument<Real>::type p,typename ScalarArgument<Real>::type pp){
  JetT<Real> c(v);for(int i=0;i<4;i++){c.d[i]=p*a.d[i];
    for(int j=0;j<4;j++)c.h[i][j]=p*a.h[i][j]+pp*a.d[i]*a.d[j];}return c;
}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> power(const JetT<Real>&a,typename ScalarArgument<Real>::type p){
  Real v=scalar_power(a.v,p);return unary(a,v,p*v/a.v,p*(p-1)*v/(a.v*a.v));
}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> inverse(const JetT<Real>&a){return power(a,-1);}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> operator/(const JetT<Real>&a,const JetT<Real>&b){return a*inverse(b);}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> sqrt(const JetT<Real>&a){return power(a,.5L);}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> exp(const JetT<Real>&a){Real v=std::exp(a.v);return unary(a,v,v,v);}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> expm1(const JetT<Real>&a){Real p=std::exp(a.v);return unary(a,std::expm1(a.v),p,p);}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> tan(const JetT<Real>&a){Real v=std::tan(a.v),p=1+v*v;return unary(a,v,p,2*v*p);}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> tanh(const JetT<Real>&a){Real v=std::tanh(a.v),p=1-v*v;return unary(a,v,p,-2*v*p);}
/* This differentiated jet is correct to FIRST order only. Used for K and
 * connection first derivatives; third spacetime derivatives are never used. */
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> diff(const JetT<Real>&a,int k){JetT<Real> c(a.d[k]);for(int i=0;i<4;i++)c.d[i]=a.h[k][i];return c;}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> determinant(const JetT<Real> g[3][3]){
 return g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[2][1])
       -g[0][1]*(g[1][0]*g[2][2]-g[1][2]*g[2][0])
       +g[0][2]*(g[1][0]*g[2][1]-g[1][1]*g[2][0]);
}
template<class Real> HISPID_GEOMETRY_INLINE bool invert_checked(const JetT<Real> g[3][3],JetT<Real> inv[3][3]){
 JetT<Real> det=determinant(g);if(!(det.v>0))return false;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  int a=(j+1)%3,b=(j+2)%3,c=(i+1)%3,d=(i+2)%3;
  inv[i][j]=(g[a][c]*g[b][d]-g[a][d]*g[b][c])/det;
 }return true;
}
inline void invert(const Jet g[3][3],Jet inv[3][3]){
 if(!invert_checked(g,inv))throw std::runtime_error("nonpositive spatial determinant");
}
template<class Real> HISPID_GEOMETRY_INLINE void connection(const JetT<Real> g[3][3],const JetT<Real> inv[3][3],JetT<Real> C[3][3][3]){
 for(int k=0;k<3;k++)for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  C[k][i][j]=0;for(int l=0;l<3;l++)C[k][i][j]=C[k][i][j]
   +JetT<Real>(.5L)*inv[k][l]*(diff(g[j][l],i+1)+diff(g[i][l],j+1)-diff(g[i][j],l+1));
 }
}
template<class Real> HISPID_GEOMETRY_INLINE Real curvature(const JetT<Real> inv[3][3],const JetT<Real> C[3][3][3]){
 Real R=0;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  Real rij=0;for(int k=0;k<3;k++){
   rij+=C[k][i][j].d[k+1]-C[k][i][k].d[j+1];
   for(int l=0;l<3;l++)rij+=C[k][i][j].v*C[l][k][l].v-C[l][i][k].v*C[k][j][l].v;
  }R+=inv[i][j].v*rij;
 }
 return R;
}
}
#endif
