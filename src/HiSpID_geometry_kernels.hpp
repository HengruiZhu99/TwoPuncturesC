#ifndef HISPID_GEOMETRY_KERNELS_HPP
#define HISPID_GEOMETRY_KERNELS_HPP
#include "HiSpID.h"
#include "HiSpID_jets.hpp"
#include "HiSpID_geometry_types.hpp"
#include "HiSpID_trumpet.hpp"
namespace hispid {
template<class Real> HISPID_GEOMETRY_INLINE int seed_geometry(const HiSpID_Hole&hole,int choice,const double *point,SeedT<Real>&s,bool stable_evaluation=false,int family=HISPID_SEED_QI){
 if(family==HISPID_SEED_TRUMPET_R0_M)return trumpet_seed_geometry(hole,choice,point,s);
 if(family!=HISPID_SEED_QI)return geometry_nonfinite;
 Real m=hole.mass,smag=point_norm<Real>(hole.spin),a=smag/m,v2=0;
 Real axis[3]={0,0,1};if(smag>0)for(int i=0;i<3;i++)axis[i]=hole.spin[i]/smag;
 for(int i=0;i<3;i++)v2+=(Real)hole.velocity[i]*hole.velocity[i];
 Real boost=1/std::sqrt(1-v2),B[4][4]={};
 B[0][0]=boost;for(int i=0;i<3;i++){
  B[0][i+1]=B[i+1][0]=-boost*hole.velocity[i];
  for(int j=0;j<3;j++)B[i+1][j+1]=(i==j?1:0)
    +(v2>0?(boost-1)*hole.velocity[i]*hole.velocity[j]/v2:0);
 }
 /* Form geometry with REST-coordinate derivative slots. The graph-slice
  * algebra below requires those derivatives; pull all jets back afterward. */
 JetT<Real> x[3];for(int i=0;i<3;i++){
  Real value=0;for(int k=0;k<3;k++)value+=B[i+1][k+1]*(point[k]-hole.center[k]);
  x[i]=JetT<Real>::variable(value,i+1);
 }
 JetT<Real> r2=0,z=0;for(int i=0;i<3;i++){r2=r2+x[i]*x[i];z=z+JetT<Real>(axis[i])*x[i];}
 if(!(r2.v>0))return geometry_puncture;
 JetT<Real> r=sqrt(r2),costh=z/r;
 JetT<Real> rb=r+JetT<Real>(m)+JetT<Real>((m*m-a*a)/4)/r;
 JetT<Real> sigma=rb*rb+JetT<Real>(a*a)*costh*costh;
 JetT<Real> sin2=JetT<Real>(1)-costh*costh;
 JetT<Real> AA=(rb*rb+JetT<Real>(a*a))*(rb*rb+JetT<Real>(a*a))
    -JetT<Real>(a*a)*(r-JetT<Real>((m*m-a*a)/4)/r)*(r-JetT<Real>((m*m-a*a)/4)/r)*sin2;
 JetT<Real> psiQI=power(sigma/r2,Real(.25L));
 JetT<Real> alpha0=(r-JetT<Real>((m*m-a*a)/4)/r)*sqrt(sigma/AA);
 JetT<Real> cross[3]={JetT<Real>(axis[1])*x[2]-JetT<Real>(axis[2])*x[1],
    JetT<Real>(axis[2])*x[0]-JetT<Real>(axis[0])*x[2],JetT<Real>(axis[0])*x[1]-JetT<Real>(axis[1])*x[0]};
 JetT<Real> beta0[3],rest[3][3],rinv[3][3];
 for(int i=0;i<3;i++){
  beta0[i]=-JetT<Real>(2*m*a)*rb/AA*cross[i];
  for(int j=0;j<3;j++)rest[i][j]=power(psiQI,4)
     *(JetT<Real>(i==j?1:0)+JetT<Real>(a*a)*(JetT<Real>(1)+JetT<Real>(2*m)*rb/sigma)/(sigma*r2)*cross[i]*cross[j]);
 }
 if(!invert_checked(rest,rinv))return geometry_determinant;
 {
 JetT<Real> grest[4][4];
 grest[0][0]=-alpha0*alpha0;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  grest[i+1][j+1]=rest[i][j];
  grest[0][i+1]=grest[0][i+1]+rest[i][j]*beta0[j];
  grest[0][0]=grest[0][0]+rest[i][j]*beta0[i]*beta0[j];
 }
 for(int i=0;i<3;i++)grest[i+1][0]=grest[0][i+1];
 // Only the spatial block is consumed. Retain the original k/l arithmetic
 // order, and reset reused output just as the former fresh g4 array did.
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  s.physical[i][j]=0;
  for(int k=0;k<4;k++)for(int l=0;l<4;l++)
   s.physical[i][j]=s.physical[i][j]+JetT<Real>(B[k][i+1]*B[l][j+1])*grest[k][l];
 }
 }
 JetT<Real> inv[3][3],C0[3][3][3];
 if(!invert_checked(s.physical,inv))return geometry_determinant;connection(rest,rinv,C0);
 /* Rest beta=omega*l, with l an axial Killing vector. Factoring dR/dr
  * and Delta analytically in (partial omega)/alpha removes the 0/0 at
  * the QI throat and the associated cancellation in K and its derivatives. */
 const JetT<Real> radial=r-JetT<Real>((m*m-a*a)/4)/r;
 JetT<Real> AAR=JetT<Real>(4)*rb*(rb*rb+JetT<Real>(a*a))-JetT<Real>(a*a)*(JetT<Real>(2)*rb-JetT<Real>(2*m))*sin2;
 JetT<Real> omegaR=-JetT<Real>(2*m*a)*(AA-rb*AAR)/(AA*AA),gradomega[3],lcov[3],K0[3][3],Kmixed[3][3];
 for(int i=0;i<3;i++){
  JetT<Real> dcosth=(JetT<Real>(axis[i])-costh*x[i]/r)/r;
  gradomega[i]=sqrt(AA/sigma)*(omegaR*x[i]/r2+JetT<Real>(4*m*a*a*a)*rb*costh*radial/(AA*AA)*dcosth);
  for(int j=0;j<3;j++)lcov[i]=lcov[i]+rest[i][j]*cross[j];
 }
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)K0[i][j]=(lcov[i]*gradomega[j]+lcov[j]*gradomega[i])/JetT<Real>(2);
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)Kmixed[i][j]=Kmixed[i][j]+K0[i][k]*rinv[k][j];
 JetT<Real> q=1,VV=0;for(int i=0;i<3;i++){
  q=q-JetT<Real>(hole.velocity[i])*beta0[i];
  for(int j=0;j<3;j++)VV=VV+JetT<Real>(hole.velocity[i]*hole.velocity[j])*rinv[i][j];
 }
 /* The boosted slice is the rest-coordinate graph t0=-v.X0. Its normal
  * is W(n0+u), tangent E_a=T_a^b e_b-alpha0*v_a*n0. Expanding its second
  * fundamental form cancels every inverse alpha0 analytically. */
 JetT<Real> T[3][3],u[3],Q[3][3],N[3],Kgraph[3][3];
 JetT<Real> denominator=q*q-alpha0*alpha0*VV;
 if(!(denominator.v>0)||!(q.v>0))return geometry_nonspacelike;
 JetT<Real> W=q/sqrt(denominator);
 for(int i=0;i<3;i++){
  for(int j=0;j<3;j++){
   T[i][j]=JetT<Real>(i==j?1:0)-JetT<Real>(hole.velocity[i])*beta0[j];
   u[i]=u[i]-alpha0/q*rinv[i][j]*JetT<Real>(hole.velocity[j]);
  }
 }
 for(int b=0;b<3;b++){
  for(int d=0;d<3;d++){
   Q[b][d]=diff(u[d],b+1);
   for(int e=0;e<3;e++)Q[b][d]=Q[b][d]-T[b][e]*Kmixed[e][d]-JetT<Real>(hole.velocity[b])*rinv[d][e]*diff(alpha0,e+1);
   for(int c=0;c<3;c++){
    Q[b][d]=Q[b][d]+alpha0*JetT<Real>(hole.velocity[b])*u[c]*Kmixed[c][d]
      -JetT<Real>(hole.velocity[b])*u[c]*diff(beta0[d],c+1);
    for(int e=0;e<3;e++)Q[b][d]=Q[b][d]+u[c]*T[b][e]*C0[d][e][c];
   }
  }
  for(int c=0;c<3;c++){
   N[b]=N[b]-JetT<Real>(hole.velocity[b])*u[c]*diff(alpha0,c+1);
   for(int e=0;e<3;e++)N[b]=N[b]-u[c]*T[b][e]*K0[e][c];
  }
 }
 for(int aidx=0;aidx<3;aidx++)for(int b=0;b<3;b++){
  Kgraph[aidx][b]=-W*alpha0*JetT<Real>(hole.velocity[aidx])*N[b];
  for(int c=0;c<3;c++)for(int d=0;d<3;d++)Kgraph[aidx][b]=Kgraph[aidx][b]-W*T[aidx][c]*rest[c][d]*Q[b][d];
 }
 s.K=0;JetT<Real> trace_roundoff;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  for(int aidx=0;aidx<3;aidx++)for(int b=0;b<3;b++)
   s.extrinsic[i][j]=s.extrinsic[i][j]+JetT<Real>(B[aidx+1][i+1]*B[b+1][j+1])*(Kgraph[aidx][b]+Kgraph[b][aidx])/JetT<Real>(2);
  if(stable_evaluation)compensated_add(s.K,trace_roundoff,inv[i][j]*s.extrinsic[i][j]);
  else s.K=s.K+inv[i][j]*s.extrinsic[i][j];
 }
 /* A stationary unboosted QI Kerr slice is exactly maximal. Retaining
  * a roundoff trace here would amplify its gradient by psi^6 at a puncture. */
 if(v2==0)s.K=0;
 s.psi=choice?power(determinant(s.physical),Real(1.0L)/12):psiQI;
 JetT<Real> normalized[3][3];
 if(stable_evaluation){
  /* Build the QI conformal spacetime metric directly. Dividing two large
   * physical-metric/psi jets creates spurious derivatives of the exact
   * constant transverse metric for a highly boosted Schwarzschild seed.
   * U is the spatial pullback including the rest shift; the lapse term
   * supplies the remaining normal part of the spacetime pullback. */
  JetT<Real> rest_normalized[3][3],U[3][3];
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   rest_normalized[i][j]=JetT<Real>(i==j?1:0)
    +JetT<Real>(a*a)*(JetT<Real>(1)+JetT<Real>(2*m)*rb/sigma)/(sigma*r2)*cross[i]*cross[j];
   U[i][j]=JetT<Real>(B[i+1][j+1])+JetT<Real>(B[0][j+1])*beta0[i];
  }
  JetT<Real> radial_normalized=r2-JetT<Real>((m*m-a*a)/4);
  JetT<Real> lapse_normalized=radial_normalized*radial_normalized/AA;
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   normalized[i][j]=-JetT<Real>(B[0][i+1]*B[0][j+1])*lapse_normalized;
   for(int k=0;k<3;k++)for(int l=0;l<3;l++)
    normalized[i][j]=normalized[i][j]+U[k][i]*rest_normalized[k][l]*U[l][j];
  }
  if(choice){
   JetT<Real> det=determinant(normalized);if(!(det.v>0))return geometry_determinant;
   JetT<Real> factor=power(det,Real(1.0L)/3);
   s.psi=psiQI*power(det,Real(1.0L)/12);
   for(int i=0;i<3;i++)for(int j=0;j<3;j++)normalized[i][j]=normalized[i][j]/factor;
  }
 }
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  s.metric[i][j]=stable_evaluation?normalized[i][j]:s.physical[i][j]/power(s.psi,4);
  s.A[i][j]=s.psi*s.psi*(s.extrinsic[i][j]-s.physical[i][j]*s.K/JetT<Real>(3));
  s.physical[i][j]=pullback_derivatives(s.physical[i][j],B);
  s.extrinsic[i][j]=pullback_derivatives(s.extrinsic[i][j],B);
  s.metric[i][j]=pullback_derivatives(s.metric[i][j],B);
  s.A[i][j]=pullback_derivatives(s.A[i][j],B);
 }
 s.psi=pullback_derivatives(s.psi,B);s.K=pullback_derivatives(s.K,B);
 return geometry_ok;
}

/* Value-only lab-frame export for both policies. K/A/mean-K Hessians are
 * intentionally absent: graph-slice differentiation provides first order.
 * Exact seeds have zero correction and no attenuation. */
template<class Real> HISPID_GEOMETRY_INLINE void seed_values(const SeedT<Real>&s,HiSpID_Point&out){
 out={};
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  int k=3*i+j;out.gamma[k]=s.physical[i][j].v;out.Kij[k]=s.extrinsic[i][j].v;
  out.conformal_metric[k]=s.metric[i][j].v;out.Atilde[k]=s.A[i][j].v;
 }
 out.psi=s.psi.v;out.mean_curvature=s.K.v;out.attenuation=1;
}
HISPID_GEOMETRY_INLINE bool finite_seed_values(const HiSpID_Point&out){
 if(!std::isfinite(out.psi)||!std::isfinite(out.mean_curvature))return false;
 for(int k=0;k<9;k++)if(!std::isfinite(out.gamma[k])||!std::isfinite(out.Kij[k])
   ||!std::isfinite(out.conformal_metric[k])||!std::isfinite(out.Atilde[k]))return false;
 return true;
}

template<class Real> HISPID_GEOMETRY_INLINE JetT<Real> inner(const JetT<Real>&r,double lo,double hi){
 if(hi<=0 || r.v>=hi)return 1;
 if(r.v<=lo)return 0;
 JetT<Real> z=(r-JetT<Real>(lo))/JetT<Real>(hi-lo);
 /* Tanh(tan) is C-infinity at both ends; saturating far past double
  * precision avoids spurious 0*infinity in its jet derivatives. */
 if(z.v<Real(.01L))return 0;if(z.v>Real(.99L))return 1;
 return (JetT<Real>(1)+tanh(tan(JetT<Real>(std::acos(Real(-1))/2)*(JetT<Real>(-1)+JetT<Real>(2)*z))))/JetT<Real>(2);
}

/* Cancel each isolated seed's singular momentum identity analytically.
 * All connections here are Levi-Civita connections of the actual metrics;
 * the modified interior correction connection is never used for the source. */
template<class Real> HISPID_GEOMETRY_INLINE int seed_sum_source(const HiSpID_Config&cfg,const SeedT<Real> s[2],
                            const JetT<Real> f[2],const JetT<Real> F[2],const BackgroundT<Real>&b,
                            JetT<Real>&trace,double out[3]){
 Real div[3]={};trace=0;
 for(int h=0;h<2;h++)if(cfg.hole[h].mass>0){
  JetT<Real> si[3][3],SC[3][3][3],dh[3][3],di[3][3],up[3][3],du[3][3];
  if(!invert_checked(s[h].metric,si))return geometry_determinant;connection(s[h].metric,si,SC);
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   dh[i][j]=(f[h]*F[h]-JetT<Real>(1))*(s[h].metric[i][j]-JetT<Real>(i==j?1:0));
   for(int other=0;other<2;other++)if(other!=h&&cfg.hole[other].mass>0)
    dh[i][j]=dh[i][j]+f[other]*F[other]*(s[other].metric[i][j]-JetT<Real>(i==j?1:0));
  }
  /* h^-1-hseed^-1 = -h^-1 (h-hseed) hseed^-1. */
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)
   for(int k=0;k<3;k++)for(int l=0;l<3;l++)di[i][j]=di[i][j]-b.inv[i][k]*dh[k][l]*si[l][j];
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   trace=trace+di[i][j]*s[h].A[i][j]; // trace_seed Aseed is identically zero
   for(int k=0;k<3;k++)for(int l=0;l<3;l++){
    up[i][j]=up[i][j]+b.inv[i][k]*b.inv[j][l]*s[h].A[k][l];
    du[i][j]=du[i][j]+(di[i][k]*b.inv[j][l]+si[i][k]*di[j][l])*s[h].A[k][l];
   }
  }
  Real Dh[3][3][3]={},DC[3][3][3]={};
  for(int d=0;d<3;d++)for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   Dh[d][i][j]=dh[i][j].d[d+1];
   for(int k=0;k<3;k++)Dh[d][i][j]-=SC[k][d][i].v*dh[k][j].v+SC[k][d][j].v*dh[i][k].v;
  }
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)
   for(int l=0;l<3;l++)DC[i][j][k]+=Real(.5L)*b.inv[i][l].v*(Dh[j][k][l]+Dh[k][j][l]-Dh[l][j][k]);
  for(int i=0;i<3;i++){
   for(int j=0;j<3;j++){
    div[i]+=Real(2.0L)/3*std::pow(s[h].psi.v,6)*si[i][j].v*s[h].K.d[j+1];
    div[i]+=du[i][j].d[j+1];
    for(int k=0;k<3;k++)div[i]+=SC[i][j][k].v*du[k][j].v+SC[j][j][k].v*du[i][k].v
      +DC[i][j][k]*up[k][j].v+DC[j][j][k]*up[i][k].v;
   }
  }
 }
 for(int i=0;i<3;i++){
  for(int j=0;j<3;j++)div[i]-=b.inv[i][j].v*trace.d[j+1]/3;
  out[i]=(double)div[i];
 }return geometry_ok;
}

template<class Real> HISPID_GEOMETRY_INLINE int background_geometry(const HiSpID_Config&cfg,const double*x,BackgroundT<Real>&b,bool stable_evaluation=false,int family=HISPID_SEED_QI){
 SeedT<Real> s[2];JetT<Real> radius[2],F[2]={1,1},f[2]={1,1},complement[2];
 b.psi=1;b.g=1;b.K=0;
 for(int h=0;h<2;h++)if(cfg.hole[h].mass>0){
  JetT<Real> r2=0;for(int i=0;i<3;i++){
   JetT<Real> dx=JetT<Real>::variable(x[i]-cfg.hole[h].center[i],i+1);r2=r2+dx*dx;
  }radius[h]=sqrt(r2);
  const int status=seed_geometry(cfg.hole[h],cfg.conformal_choice,x,s[h],stable_evaluation,family);if(status)return status;
  if(cfg.far_radius>0){JetT<Real> exponent=-power(radius[h]/JetT<Real>(cfg.far_radius),4);
   F[h]=exp(exponent);complement[h]=stable_evaluation?-expm1(exponent):JetT<Real>(1)-F[h];}
  b.psi=b.psi+F[h]*(s[h].psi-JetT<Real>(1));
  b.far_correction=b.far_correction+complement[h]*(s[h].psi-JetT<Real>(1));
  b.g=b.g*inner(radius[h],cfg.inner_min[h],cfg.inner_max[h]);
 }
 for(int h=0;h<2;h++)if(cfg.omega[h]>0&&cfg.hole[1-h].mass>0){
  JetT<Real> exponent=-power(radius[1-h]/JetT<Real>(cfg.omega[h]),cfg.attenuation_power);
  f[h]=stable_evaluation?-expm1(exponent):JetT<Real>(1)-exp(exponent);
 }
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  b.metric[i][j]=JetT<Real>(i==j?1:0);b.M[i][j]=0;
  for(int h=0;h<2;h++)if(cfg.hole[h].mass>0){
   b.metric[i][j]=b.metric[i][j]+f[h]*F[h]*(s[h].metric[i][j]-JetT<Real>(i==j?1:0));
   b.M[i][j]=b.M[i][j]+s[h].A[i][j];
  }
 }
 for(int h=0;h<2;h++)if(cfg.hole[h].mass>0)b.K=b.K+f[h]*F[h]*s[h].K;
 if(!invert_checked(b.metric,b.inv))return geometry_determinant;connection(b.metric,b.inv,b.C);
 JetT<Real> trace;const int status=seed_sum_source(cfg,s,f,F,b,trace,b.divM);if(status)return status;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)b.M[i][j]=b.M[i][j]-b.metric[i][j]*trace/JetT<Real>(3);
 b.R=(double)curvature(b.inv,b.C);
 b.lapPsi=geometry_laplacian(b.inv,b.C,b.psi);
 if(cfg.inner_flatten){
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)
   b.opmetric[i][j]=JetT<Real>(i==j?1:0)+b.g*(b.metric[i][j]-JetT<Real>(i==j?1:0));
  if(!invert_checked(b.opmetric,b.opinv))return geometry_determinant;
  for(int k=0;k<3;k++)for(int i=0;i<3;i++)for(int j=0;j<3;j++)b.opC[k][i][j]=b.g*b.C[k][i][j];
 }else{
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){b.opmetric[i][j]=b.metric[i][j];b.opinv[i][j]=b.inv[i][j];
   for(int k=0;k<3;k++)b.opC[k][i][j]=b.C[k][i][j];}
 }
 if(!std::isfinite(b.R)||!std::isfinite(b.lapPsi))return geometry_nonfinite;
 return geometry_ok;
}

template<class Real> HISPID_GEOMETRY_INLINE double geometry_laplacian(const JetT<Real> inv[3][3],const JetT<Real> C[3][3][3],const JetT<Real>&u){
 Real out=0;for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  Real h=u.h[i+1][j+1];for(int k=0;k<3;k++)h-=C[k][i][j].v*u.d[k+1];
  out+=inv[i][j].v*h;
 }
 return (double)out;
}

}
#endif
