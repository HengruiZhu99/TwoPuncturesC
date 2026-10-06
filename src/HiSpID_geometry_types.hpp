#ifndef HISPID_GEOMETRY_TYPES_HPP
#define HISPID_GEOMETRY_TYPES_HPP
#include "HiSpID.h"
#include "HiSpID_jets.hpp"
namespace hispid {
enum GeometryStatus {geometry_ok=0,geometry_puncture=1,geometry_determinant=2,geometry_nonspacelike=3,geometry_nonfinite=4};
template<class Real> struct SeedT {JetT<Real> metric[3][3],physical[3][3],extrinsic[3][3],A[3][3],psi,K;};
template<class Real> struct BackgroundT {
 JetT<Real> metric[3][3],inv[3][3],C[3][3][3],M[3][3],psi,K,g,far_correction;
 JetT<Real> opmetric[3][3],opinv[3][3],opC[3][3][3];
 double R,lapPsi,divM[3];
};
using Seed=SeedT<long double>;
using Background=BackgroundT<long double>;
template<class Real> HISPID_GEOMETRY_INLINE double geometry_laplacian(const JetT<Real> inv[3][3],const JetT<Real> C[3][3][3],const JetT<Real>&u);
template<class Real> HISPID_GEOMETRY_INLINE Real point_norm(const double *a){return std::sqrt((Real)a[0]*a[0]+(Real)a[1]*a[1]+(Real)a[2]*a[2]);}
template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> pullback_derivatives(const JetT<Real>&f,const Real B[4][4]){
 JetT<Real> out(f.v);
 for(int a=0;a<4;a++)for(int k=0;k<4;k++)out.d[a]+=f.d[k]*B[k][a];
 for(int a=0;a<4;a++)for(int b=0;b<4;b++)for(int k=0;k<4;k++)for(int l=0;l<4;l++)
  out.h[a][b]+=f.h[k][l]*B[k][a]*B[l][b];
 return out;
}

}
#endif
