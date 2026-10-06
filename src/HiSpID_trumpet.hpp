#ifndef HISPID_TRUMPET_HPP
#define HISPID_TRUMPET_HPP
#include "HiSpID_geometry_kernels.hpp"

namespace hispid {
// R0=M stationary Kerr trumpet, DBM arXiv:1409.1887, Eq.12.
// Cartesian algebra adapted from AthenaK project/kerr-trumpet-spin09,
// commit 81c7d803, src/coordinates/kerr_trumpet.hpp. Spin is a rest vector.
// Full spacetime pullback; metric/scalar jets are second order, K/A first
// order only (as in the existing SeedT contract). No K Hessian is supplied.
template<class Real> HISPID_GEOMETRY_INLINE int trumpet_seed_geometry(
    const HiSpID_Hole &hole, int choice, const double *point, SeedT<Real> &s,
    Real *slice_margin=nullptr) {
  using J=JetT<Real>;
  s=SeedT<Real>();
  const Real m=hole.mass, smag=point_norm<Real>(hole.spin);
  Real v2=0;
  for(int i=0;i<3;i++) v2+=Real(hole.velocity[i])*hole.velocity[i];
  if(!(m>0) || !(smag<m*m) || !(v2<1) || (choice!=0 && choice!=1))
    return geometry_nonfinite;
  const Real a=smag/m, gap=std::sqrt(m*m-a*a), boost=1/std::sqrt(1-v2);
  Real axis[3]={0,0,1}, B[4][4]={};
  if(smag>0) for(int i=0;i<3;i++) axis[i]=hole.spin[i]/smag;
  B[0][0]=boost;
  for(int i=0;i<3;i++) {
    B[0][i+1]=B[i+1][0]=-boost*hole.velocity[i];
    for(int j=0;j<3;j++) B[i+1][j+1]=Real(i==j)
      +(v2>0 ? (boost-1)*hole.velocity[i]*hole.velocity[j]/v2 : 0);
  }
  // Coordinates carry LAB spacetime derivatives, including moving-center dt.
  J x[3], r2, z;
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) x[i].v+=B[i+1][j+1]*(point[j]-hole.center[j]);
    for(int mu=0;mu<4;mu++) x[i].d[mu]=B[i+1][mu];
    r2=r2+x[i]*x[i]; z=z+J(axis[i])*x[i];
  }
  if(!(r2.v>0)) return geometry_puncture;
  const J r=sqrt(r2), R=r+J(m), ct=z/r;
  const J Sigma=R*R+J(a*a)*ct*ct, A=R*R+J(a*a);
  const J X=A*A-J(a*a)*(r2-z*z);
  const J alpha0=r*sqrt(Sigma/X);
  J cross[3]={J(axis[1])*x[2]-J(axis[2])*x[1],
              J(axis[2])*x[0]-J(axis[0])*x[2],
              J(axis[0])*x[1]-J(axis[1])*x[0]};
  J rest[3][3], ri[3][3], beta0[3], g0[4][4];
  for(int i=0;i<3;i++) {
    beta0[i]=(J(gap)*A*x[i]-J(a)*(J(2*m)*r+J(m*m+a*a))*cross[i])/X;
    for(int j=0;j<3;j++) rest[i][j]=(Sigma*J(i==j)
      +J(a*a)*(J(1)+J(2*m)*R/Sigma)*cross[i]*cross[j]/r2
      -J(a*gap)*(x[i]*cross[j]+cross[i]*x[j])/r2)/r2;
  }
  if(!invert_checked(rest,ri)) return geometry_determinant;
  J q(1), vv;
  for(int i=0;i<3;i++) {
    q=q-J(hole.velocity[i])*beta0[i];
    for(int j=0;j<3;j++) vv=vv+J(Real(hole.velocity[i])*hole.velocity[j])*ri[i][j];
  }
  const J margin=q*q-alpha0*alpha0*vv;
  if(slice_margin) *slice_margin=margin.v;
  if(!(margin.v>0) || !(q.v>0)) return geometry_nonspacelike;
  g0[0][0]=-alpha0*alpha0;
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) {
    g0[i+1][j+1]=rest[i][j];
    g0[0][i+1]=g0[0][i+1]+rest[i][j]*beta0[j];
    g0[0][0]=g0[0][0]+rest[i][j]*beta0[i]*beta0[j];
  }
  for(int i=0;i<3;i++) g0[i+1][0]=g0[0][i+1];
  J g[4][4];
  for(int mu=0;mu<4;mu++) for(int nu=0;nu<4;nu++)
    for(int k=0;k<4;k++) for(int l=0;l<4;l++)
      g[mu][nu]=g[mu][nu]+J(B[k][mu]*B[l][nu])*g0[k][l];
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) s.physical[i][j]=g[i+1][j+1];
  // Sylvester criteria, not determinant alone, enforce a spacelike slice.
  if(!(g[1][1].v>0) || !(g[1][1].v*g[2][2].v-g[1][2].v*g[1][2].v>0))
    return geometry_nonspacelike;
  J inv[3][3], beta[3];
  if(!invert_checked(s.physical,inv)) return geometry_nonspacelike;
  const J alpha=alpha0/(J(boost)*sqrt(margin));
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) beta[i]=beta[i]+inv[i][j]*g[0][j+1];
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) {
    J lie=-diff(s.physical[i][j],0);
    for(int k=0;k<3;k++) lie=lie+beta[k]*diff(s.physical[i][j],k+1)
      +s.physical[k][j]*diff(beta[k],i+1)+s.physical[i][k]*diff(beta[k],j+1);
    s.extrinsic[i][j]=lie/(J(2)*alpha);
    s.K=s.K+inv[i][j]*s.extrinsic[i][j];
  }
  s.psi=choice ? power(determinant(s.physical),Real(1)/12) : power(Sigma/r2,Real(.25));
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) {
    s.metric[i][j]=s.physical[i][j]/power(s.psi,4);
    // CTT covariant A = psi^2 (K-gamma K/3), not Z4c's psi^-4 TF(K).
    s.A[i][j]=s.psi*s.psi*(s.extrinsic[i][j]-s.physical[i][j]*s.K/J(3));
    if(!std::isfinite(s.physical[i][j].v) || !std::isfinite(s.extrinsic[i][j].v))
      return geometry_nonfinite;
  }
  return std::isfinite(s.psi.v) && s.psi.v>0 ? geometry_ok : geometry_nonfinite;
}
}
#endif
