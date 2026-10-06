#ifndef HISPID_SYMMETRY_HPP
#define HISPID_SYMMETRY_HPP
#include <cmath>
#ifdef PUNCTURES_KOKKOS
#include <Kokkos_Core.hpp>
#define HISPID_SYMMETRY_INLINE KOKKOS_INLINE_FUNCTION
#else
#define HISPID_SYMMETRY_INLINE inline
#endif
namespace hispid {
// Orthogonal no-swirl SO(2) projectors about the local puncture x axis.
// Unknowns use orthonormal Fourier slots; equations use physical phi nodes.
HISPID_SYMMETRY_INLINE void axisymmetric_unknown_row(double *v,int row,int stride,int np) {
 const double u=v[4*row],bx=v[4*row+1];
 const double br=.5*(v[4*(row+stride)+2]+v[4*(row+stride*(np/2+1))+3]);
 for(int k=0;k<np;k++)for(int c=0;c<4;c++)v[4*(row+stride*k)+c]=0;
 v[4*row]=u;v[4*row+1]=bx;
 v[4*(row+stride)+2]=br;v[4*(row+stride*(np/2+1))+3]=br;
}
HISPID_SYMMETRY_INLINE void axisymmetric_equation_row(double *v,int row,int stride,int np) {
 constexpr double tau=6.283185307179586476925286766559;
 double u=0,bx=0,br=0;
 for(int k=0;k<np;k++){
  const int p=4*(row+stride*k);const double phi=tau*k/np;
  u+=v[p]/np;bx+=v[p+1]/np;
  br+=(std::cos(phi)*v[p+2]+std::sin(phi)*v[p+3])/np;
 }
 for(int k=0;k<np;k++){
  const int p=4*(row+stride*k);const double phi=tau*k/np;
  v[p]=u;v[p+1]=bx;v[p+2]=br*std::cos(phi);v[p+3]=br*std::sin(phi);
 }
}
}
#undef HISPID_SYMMETRY_INLINE
#endif
