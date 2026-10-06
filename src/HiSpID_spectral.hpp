#ifndef HISPID_SPECTRAL_HPP
#define HISPID_SPECTRAL_HPP
#include <array>
#include <vector>
#include <cmath>
#include "HiSpID_jets.hpp"
#ifndef HISPID_ANALYTIC_MATRICES
#define HISPID_ANALYTIC_MATRICES 0
#endif
extern "C" {
#include "TwoPunctures.h"
}
namespace hispid {
// Chebyshev-root barycentric entries, using a trigonometric node difference
// to avoid subtracting almost equal cosines near the endpoints. The second
// derivative identity is D2_ij=2 D_ij (D_ii-1/(x_i-x_j)), i!=j.
// Host builders use long double; device cache builders use double.
template<class Real> HISPID_GEOMETRY_INLINE Real chebyshev_offdiagonal(int N,int i,int j,int order){
 const Real pi=std::acos(Real(-1)),ti=pi*(Real(i)+Real(.5))/Real(N),tj=pi*(Real(j)+Real(.5))/Real(N);
 const Real si=std::sin(ti),dx=Real(2)*std::sin((ti+tj)/Real(2))*std::sin((ti-tj)/Real(2));
 const Real d=(std::abs(i-j)%2?Real(-1):Real(1))*std::sin(tj)/(si*dx);
 return order==1?d:Real(2)*d*(-std::cos(ti)/(Real(2)*si*si)-Real(1)/dx);
}

/* Differentiation matrices for exactly the inherited interpolant. Reusing
 * them avoids repeated slow trigonometric transforms in every Krylov call.
 * Fourier D2 is formed separately: D*D would lose the cosine Nyquist mode. */
struct SpectralDerivatives {
 std::array<int,3> n{};
 std::array<std::vector<double>,3> D,D2;
 void initialize(const int*shape){
  const double pi=std::acos(-1.0);
  for(int axis=0;axis<3;axis++){
   const int N=n[axis]=shape[axis];D[axis].assign(N*N,0);D2[axis].assign(N*N,0);
   if(axis<2){
#if HISPID_ANALYTIC_MATRICES
    for(int i=0;i<N;i++){
     long double diagonal=0,diagonal2=0;
     for(int j=0;j<N;j++)if(j!=i){
      const long double d=chebyshev_offdiagonal<long double>(N,i,j,1),d2=chebyshev_offdiagonal<long double>(N,i,j,2);
      D[axis][i*N+j]=double(d);D2[axis][i*N+j]=double(d2);diagonal-=d;diagonal2-=d2;
     }
     D[axis][i*N+i]=double(diagonal);D2[axis][i*N+i]=double(diagonal2);
    }
#else
    std::vector<double>x(N),w(N);
    for(int i=0;i<N;i++){const double t=pi*(i+.5)/N;x[i]=-std::cos(t);w[i]=(i%2?-1:1)*std::sin(t);}
    for(int i=0;i<N;i++)for(int j=0;j<N;j++)if(i!=j){
     D[axis][i*N+j]=w[j]/(w[i]*(x[i]-x[j]));D[axis][i*N+i]-=D[axis][i*N+j];
    }
    for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int k=0;k<N;k++)
     D2[axis][i*N+j]+=D[axis][i*N+k]*D[axis][k*N+j];
    // The interpolant differentiates constants to zero. Preserve that
    // identity when forming the rounded second-derivative matrix as well.
    for(int i=0;i<N;i++){
     double diagonal=0;for(int j=0;j<N;j++)if(i!=j)diagonal-=D2[axis][i*N+j];
     D2[axis][i*N+i]=diagonal;
    }
#endif
   }else{
    for(int i=0;i<N;i++){
     D2[axis][i*N+i]=-(N*N+2.0)/12;
     for(int j=0;j<N;j++)if(i!=j){
      const double t=pi*(i-j)/N,sign=(std::abs(i-j)%2?-1:1);
      D[axis][i*N+j]=.5*sign/std::tan(t);
      D2[axis][i*N+j]=-.5*sign/(std::sin(t)*std::sin(t));
     }
    }
   }
  }
 }
 void differentiate(int axis,const std::vector<double>&matrix,int nv,const double*in,double*out)const{
  const int stride=axis==0?1:axis==1?n[0]:n[0]*n[1],N=n[axis],points=n[0]*n[1]*n[2];
  // Sum off-diagonal differences instead of large nearly cancelling terms.
  // This is the same polynomial derivative and annihilates constants exactly.
  for(int block=0;block<points/(N*stride);block++)for(int offset=0;offset<stride;offset++)for(int i=0;i<N;i++){
   const int base=block*N*stride+offset,p=base+i*stride;
   for(int v=0;v<nv;v++)out[nv*p+v]=0;
   for(int j=0;j<N;j++)if(j!=i){
    const double c=matrix[i*N+j];const double*row=in+nv*(base+j*stride);
    for(int v=0;v<nv;v++)out[nv*p+v]+=c*(row[v]-in[nv*p+v]);
   }
  }
 }
 void along(int axis,const std::vector<double>&matrix,int nv,const double*in,double*out)const{
  const int stride=axis==0?1:axis==1?n[0]:n[0]*n[1],N=n[axis],points=n[0]*n[1]*n[2];
  // Reuse each meridional/azimuthal input line while it is in cache.
  // The column summation order for every output is unchanged.
  for(int block=0;block<points/(N*stride);block++)for(int offset=0;offset<stride;offset++)for(int i=0;i<N;i++){
   const int base=block*N*stride+offset,p=base+i*stride;
   for(int v=0;v<nv;v++)out[nv*p+v]=0;
   for(int j=0;j<N;j++){
    const double c=matrix[i*N+j];const double*row=in+nv*(base+j*stride);
    for(int v=0;v<nv;v++)out[nv*p+v]+=c*row[v];
   }
  }
 }
 void apply(int nv,derivs*w)const{
  along(0,D[0],nv,w->d0,w->d1);along(1,D[1],nv,w->d0,w->d2);along(2,D[2],nv,w->d0,w->d3);
  along(0,D2[0],nv,w->d0,w->d11);along(1,D2[1],nv,w->d0,w->d22);along(2,D2[2],nv,w->d0,w->d33);
  along(0,D[0],nv,w->d2,w->d12);along(0,D[0],nv,w->d3,w->d13);along(1,D[1],nv,w->d3,w->d23);
 }
};
}
#endif
