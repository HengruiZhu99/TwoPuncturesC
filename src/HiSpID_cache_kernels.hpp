#ifndef HISPID_CACHE_KERNELS_HPP
#define HISPID_CACHE_KERNELS_HPP
#include "HiSpID_geometry_kernels.hpp"
#include <array>
namespace hispid {
#ifdef PUNCTURES_KOKKOS
using CachedFields=Kokkos::Array<Kokkos::Array<double,10>,4>;
#else
using CachedFields=std::array<std::array<double,10>,4>;
#endif
HISPID_GEOMETRY_INLINE int hessian_i(int d){const int values[6]={1,1,1,2,2,3};return values[d];}
HISPID_GEOMETRY_INLINE int hessian_j(int d){const int values[6]={1,2,3,2,3,3};return values[d];}
struct Cached {
 double psi,R,K,g,lapPsi,divM[3],gradK[3],M[9],inv[9];
 double lap[10],vec[3][3][10],L[9][3][4];
 double far_correction[10];
 double weight;
};
/* First-order tensor algebra for divergence of L. The scalar/vector fields
 * carry second derivatives, but background metric/connection need only
 * first derivatives here. This is also the explicit nonmetric-compatible
 * definition used for Eq27/28 inside attenuation zones. */
template<class Real> HISPID_GEOMETRY_INLINE void L_and_div(const JetT<Real> metric[3][3],const JetT<Real> inv[3][3],const JetT<Real> C[3][3][3],
               const CachedFields&u,double L[3][3],double*out){
 double db[3][3]={},ddb[3][3][3]={},div=0,ddiv[3]={};
 double bu[3][3][3]={};
 for(int k=0;k<3;k++)for(int d=0;d<6;d++)bu[k][hessian_i(d)-1][hessian_j(d)-1]=bu[k][hessian_j(d)-1][hessian_i(d)-1]=u[k+1][4+d];
 for(int i=0;i<3;i++)for(int k=0;k<3;k++){
  db[i][k]=u[k+1][i+1];for(int l=0;l<3;l++)db[i][k]+=(double)C[k][i][l].v*u[l+1][0];
  for(int d=0;d<3;d++){
   ddb[d][i][k]=bu[k][i][d];for(int l=0;l<3;l++)ddb[d][i][k]+=(double)C[k][i][l].d[d+1]*u[l+1][0]+(double)C[k][i][l].v*u[l+1][d+1];
  }
  if(i==k){div+=db[i][k];for(int d=0;d<3;d++)ddiv[d]+=ddb[d][i][k];}
 }
 double dL[3][3][3]={};
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  L[i][j]=-2.0/3*(double)metric[i][j].v*div;
  for(int k=0;k<3;k++)L[i][j]+=(double)metric[i][k].v*db[j][k]+(double)metric[j][k].v*db[i][k];
  for(int d=0;d<3;d++){
   dL[d][i][j]=-2.0/3*((double)metric[i][j].d[d+1]*div+(double)metric[i][j].v*ddiv[d]);
   for(int k=0;k<3;k++)dL[d][i][j]+=(double)metric[i][k].d[d+1]*db[j][k]+(double)metric[i][k].v*ddb[d][j][k]
     +(double)metric[j][k].d[d+1]*db[i][k]+(double)metric[j][k].v*ddb[d][i][k];
  }
 }
 if(!out)return;
 double up[3][3]={},dup[3][3][3]={};
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)for(int l=0;l<3;l++){
  up[i][j]+=(double)(inv[i][k].v*inv[j][l].v)*L[k][l];
  for(int d=0;d<3;d++)dup[d][i][j]+=(double)(inv[i][k].v*inv[j][l].v)*dL[d][k][l]
    +(double)(inv[i][k].d[d+1]*inv[j][l].v+inv[i][k].v*inv[j][l].d[d+1])*L[k][l];
 }
 for(int i=0;i<3;i++){
  out[i]=0;for(int j=0;j<3;j++){
   out[i]+=dup[j][i][j];for(int k=0;k<3;k++)out[i]+=(double)C[i][j][k].v*up[k][j]+(double)C[j][j][k].v*up[i][k];
  }
 }
}
template<class Real> HISPID_GEOMETRY_INLINE void cache(const BackgroundT<Real>&b,Cached&c){
 c=Cached{};c.psi=b.psi.v;c.R=b.R;c.K=b.K.v;c.g=b.g.v;c.lapPsi=b.lapPsi;
 c.far_correction[0]=b.far_correction.v;
 for(int d=1;d<4;d++)c.far_correction[d]=b.far_correction.d[d];
 for(int d=0;d<6;d++)c.far_correction[d+4]=b.far_correction.h[hessian_i(d)][hessian_j(d)];
 for(int i=0;i<3;i++){
  c.divM[i]=b.divM[i];for(int j=0;j<3;j++){
   c.gradK[i]+=(double)(b.inv[i][j].v*b.K.d[j+1]);
   c.M[3*i+j]=b.M[i][j].v;c.inv[3*i+j]=b.inv[i][j].v;
  }
 }
 for(int d=0;d<10;d++){
  JetT<Real> u;if(d==0)u.v=1;else if(d<4)u.d[d]=1;else u.h[hessian_i(d-4)][hessian_j(d-4)]=u.h[hessian_j(d-4)][hessian_i(d-4)]=1;
  c.lap[d]=geometry_laplacian(b.opinv,b.opC,u);
  for(int k=0;k<3;k++){
   CachedFields f{};f[k+1][d]=1;double L[3][3],div[3];
   L_and_div(b.opmetric,b.opinv,b.opC,f,L,div);
   for(int i=0;i<3;i++)c.vec[i][k][d]=div[i];
   if(d<4){L_and_div(b.metric,b.inv,b.C,f,L,nullptr);for(int i=0;i<3;i++)for(int j=0;j<3;j++)c.L[3*i+j][k][d]=L[i][j];}
  }
 }
}

HISPID_GEOMETRY_INLINE bool finite_cache(const Cached&c){
 if(!std::isfinite(c.psi)||!std::isfinite(c.R)||!std::isfinite(c.K)||!std::isfinite(c.g)||!std::isfinite(c.lapPsi)||!std::isfinite(c.weight))return false;
 for(int i=0;i<3;i++)if(!std::isfinite(c.divM[i])||!std::isfinite(c.gradK[i]))return false;
 for(int i=0;i<9;i++)if(!std::isfinite(c.M[i])||!std::isfinite(c.inv[i]))return false;
 for(int d=0;d<10;d++){
  if(!std::isfinite(c.lap[d])||!std::isfinite(c.far_correction[d]))return false;
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)if(!std::isfinite(c.vec[i][j][d]))return false;
 }
 for(int i=0;i<9;i++)for(int j=0;j<3;j++)for(int d=0;d<4;d++)if(!std::isfinite(c.L[i][j][d]))return false;
 return true;
}
}
#endif
