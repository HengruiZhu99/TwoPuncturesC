#include "HiSpID_internal.hpp"
#include <algorithm>
#include <cstring>
#include <limits>
namespace hispid {
thread_local std::string last_error;
static long double norm(const double *a){return std::sqrt((long double)a[0]*a[0]+(long double)a[1]*a[1]+(long double)a[2]*a[2]);}
static Jet pullback_derivatives(const Jet&f,const long double B[4][4]){
 Jet out(f.v);
 for(int a=0;a<4;a++)for(int k=0;k<4;k++)out.d[a]+=f.d[k]*B[k][a];
 for(int a=0;a<4;a++)for(int b=0;b<4;b++)for(int k=0;k<4;k++)for(int l=0;l<4;l++)
  out.h[a][b]+=f.h[k][l]*B[k][a]*B[l][b];
 return out;
}

void seed(const HiSpID_Hole&hole,int choice,const double *point,Seed&s){
 long double m=hole.mass,smag=norm(hole.spin),a=smag/m,v2=0;
 long double axis[3]={0,0,1};if(smag>0)for(int i=0;i<3;i++)axis[i]=hole.spin[i]/smag;
 for(int i=0;i<3;i++)v2+=(long double)hole.velocity[i]*hole.velocity[i];
 long double boost=1/std::sqrt(1-v2),B[4][4]={};
 B[0][0]=boost;for(int i=0;i<3;i++){
  B[0][i+1]=B[i+1][0]=-boost*hole.velocity[i];
  for(int j=0;j<3;j++)B[i+1][j+1]=(i==j?1:0)
    +(v2>0?(boost-1)*hole.velocity[i]*hole.velocity[j]/v2:0);
 }
 /* Form geometry with REST-coordinate derivative slots. The graph-slice
  * algebra below requires those derivatives; pull all jets back afterward. */
 Jet x[3];for(int i=0;i<3;i++){
  long double value=0;for(int k=0;k<3;k++)value+=B[i+1][k+1]*(point[k]-hole.center[k]);
  x[i]=Jet::variable(value,i+1);
 }
 Jet r2=0,z=0;for(int i=0;i<3;i++){r2=r2+x[i]*x[i];z=z+Jet(axis[i])*x[i];}
 if(!(r2.v>0))throw std::runtime_error("sampling at unsmoothed puncture");
 Jet r=sqrt(r2),costh=z/r;
 Jet rb=r+Jet(m)+Jet((m*m-a*a)/4)/r;
 Jet sigma=rb*rb+Jet(a*a)*costh*costh;
 Jet sin2=Jet(1)-costh*costh;
 Jet AA=(rb*rb+Jet(a*a))*(rb*rb+Jet(a*a))
    -Jet(a*a)*(r-Jet((m*m-a*a)/4)/r)*(r-Jet((m*m-a*a)/4)/r)*sin2;
 Jet psiQI=power(sigma/r2,.25L);
 Jet alpha0=(r-Jet((m*m-a*a)/4)/r)*sqrt(sigma/AA);
 Jet cross[3]={Jet(axis[1])*x[2]-Jet(axis[2])*x[1],
    Jet(axis[2])*x[0]-Jet(axis[0])*x[2],Jet(axis[0])*x[1]-Jet(axis[1])*x[0]};
 Jet beta0[3],rest[3][3],rinv[3][3],grest[4][4],g4[4][4];
 for(int i=0;i<3;i++){
  beta0[i]=-Jet(2*m*a)*rb/AA*cross[i];
  for(int j=0;j<3;j++)rest[i][j]=power(psiQI,4)
     *(Jet(i==j?1:0)+Jet(a*a)*(Jet(1)+Jet(2*m)*rb/sigma)/(sigma*r2)*cross[i]*cross[j]);
 }
 invert(rest,rinv);
 grest[0][0]=-alpha0*alpha0;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  grest[i+1][j+1]=rest[i][j];
  grest[0][i+1]=grest[0][i+1]+rest[i][j]*beta0[j];
  grest[0][0]=grest[0][0]+rest[i][j]*beta0[i]*beta0[j];
 }
 for(int i=0;i<3;i++)grest[i+1][0]=grest[0][i+1];
 for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int k=0;k<4;k++)for(int l=0;l<4;l++)
  g4[i][j]=g4[i][j]+Jet(B[k][i]*B[l][j])*grest[k][l];
 Jet inv[3][3],C0[3][3][3];
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)s.physical[i][j]=g4[i+1][j+1];
 invert(s.physical,inv);connection(rest,rinv,C0);
 /* Rest beta=omega*l, with l an axial Killing vector. Factoring dR/dr
  * and Delta analytically in (partial omega)/alpha removes the 0/0 at
  * the QI throat and the associated cancellation in K and its derivatives. */
 const Jet radial=r-Jet((m*m-a*a)/4)/r;
 Jet AAR=Jet(4)*rb*(rb*rb+Jet(a*a))-Jet(a*a)*(Jet(2)*rb-Jet(2*m))*sin2;
 Jet omegaR=-Jet(2*m*a)*(AA-rb*AAR)/(AA*AA),gradomega[3],lcov[3],K0[3][3],Kmixed[3][3];
 for(int i=0;i<3;i++){
  Jet dcosth=(Jet(axis[i])-costh*x[i]/r)/r;
  gradomega[i]=sqrt(AA/sigma)*(omegaR*x[i]/r2+Jet(4*m*a*a*a)*rb*costh*radial/(AA*AA)*dcosth);
  for(int j=0;j<3;j++)lcov[i]=lcov[i]+rest[i][j]*cross[j];
 }
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)K0[i][j]=(lcov[i]*gradomega[j]+lcov[j]*gradomega[i])/Jet(2);
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)Kmixed[i][j]=Kmixed[i][j]+K0[i][k]*rinv[k][j];
 Jet q=1,VV=0;for(int i=0;i<3;i++){
  q=q-Jet(hole.velocity[i])*beta0[i];
  for(int j=0;j<3;j++)VV=VV+Jet(hole.velocity[i]*hole.velocity[j])*rinv[i][j];
 }
 /* The boosted slice is the rest-coordinate graph t0=-v.X0. Its normal
  * is W(n0+u), tangent E_a=T_a^b e_b-alpha0*v_a*n0. Expanding its second
  * fundamental form cancels every inverse alpha0 analytically. */
 Jet T[3][3],u[3],Q[3][3],N[3],Kgraph[3][3];
 Jet denominator=q*q-alpha0*alpha0*VV;
 if(!(denominator.v>0)||!(q.v>0))throw std::runtime_error("boosted slice is not spacelike");
 Jet W=q/sqrt(denominator);
 for(int i=0;i<3;i++){
  for(int j=0;j<3;j++){
   T[i][j]=Jet(i==j?1:0)-Jet(hole.velocity[i])*beta0[j];
   u[i]=u[i]-alpha0/q*rinv[i][j]*Jet(hole.velocity[j]);
  }
 }
 for(int b=0;b<3;b++){
  for(int d=0;d<3;d++){
   Q[b][d]=diff(u[d],b+1);
   for(int e=0;e<3;e++)Q[b][d]=Q[b][d]-T[b][e]*Kmixed[e][d]-Jet(hole.velocity[b])*rinv[d][e]*diff(alpha0,e+1);
   for(int c=0;c<3;c++){
    Q[b][d]=Q[b][d]+alpha0*Jet(hole.velocity[b])*u[c]*Kmixed[c][d]
      -Jet(hole.velocity[b])*u[c]*diff(beta0[d],c+1);
    for(int e=0;e<3;e++)Q[b][d]=Q[b][d]+u[c]*T[b][e]*C0[d][e][c];
   }
  }
  for(int c=0;c<3;c++){
   N[b]=N[b]-Jet(hole.velocity[b])*u[c]*diff(alpha0,c+1);
   for(int e=0;e<3;e++)N[b]=N[b]-u[c]*T[b][e]*K0[e][c];
  }
 }
 for(int aidx=0;aidx<3;aidx++)for(int b=0;b<3;b++){
  Kgraph[aidx][b]=-W*alpha0*Jet(hole.velocity[aidx])*N[b];
  for(int c=0;c<3;c++)for(int d=0;d<3;d++)Kgraph[aidx][b]=Kgraph[aidx][b]-W*T[aidx][c]*rest[c][d]*Q[b][d];
 }
 s.K=0;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  for(int aidx=0;aidx<3;aidx++)for(int b=0;b<3;b++)
   s.extrinsic[i][j]=s.extrinsic[i][j]+Jet(B[aidx+1][i+1]*B[b+1][j+1])*(Kgraph[aidx][b]+Kgraph[b][aidx])/Jet(2);
  s.K=s.K+inv[i][j]*s.extrinsic[i][j];
 }
 /* A stationary unboosted QI Kerr slice is exactly maximal. Retaining
  * a roundoff trace here would amplify its gradient by psi^6 at a puncture. */
 if(v2==0)s.K=0;
 s.psi=choice?power(determinant(s.physical),1.0L/12):psiQI;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  s.metric[i][j]=s.physical[i][j]/power(s.psi,4);
  s.A[i][j]=s.psi*s.psi*(s.extrinsic[i][j]-s.physical[i][j]*s.K/Jet(3));
  s.physical[i][j]=pullback_derivatives(s.physical[i][j],B);
  s.extrinsic[i][j]=pullback_derivatives(s.extrinsic[i][j],B);
  s.metric[i][j]=pullback_derivatives(s.metric[i][j],B);
  s.A[i][j]=pullback_derivatives(s.A[i][j],B);
 }
 s.psi=pullback_derivatives(s.psi,B);s.K=pullback_derivatives(s.K,B);
}

static Jet inner(const Jet&r,double lo,double hi){
 if(hi<=0 || r.v>=hi)return 1;
 if(r.v<=lo)return 0;
 Jet z=(r-Jet(lo))/Jet(hi-lo);
 /* Tanh(tan) is C-infinity at both ends; saturating far past double
  * precision avoids spurious 0*infinity in its jet derivatives. */
 if(z.v<.01L)return 0;if(z.v>.99L)return 1;
 return (Jet(1)+tanh(tan(Jet(acosl(-1.0L)/2)*(Jet(-1)+Jet(2)*z))))/Jet(2);
}

/* Cancel each isolated seed's singular momentum identity analytically.
 * All connections here are Levi-Civita connections of the actual metrics;
 * the modified interior correction connection is never used for the source. */
static void seed_sum_source(const HiSpID_Config&cfg,const Seed s[2],
                            const Jet f[2],const Jet F[2],const Background&b,
                            Jet&trace,double out[3]){
 long double div[3]={};trace=0;
 for(int h=0;h<2;h++)if(cfg.hole[h].mass>0){
  Jet si[3][3],SC[3][3][3],dh[3][3],di[3][3],up[3][3],du[3][3];
  invert(s[h].metric,si);connection(s[h].metric,si,SC);
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   dh[i][j]=(f[h]*F[h]-Jet(1))*(s[h].metric[i][j]-Jet(i==j?1:0));
   for(int other=0;other<2;other++)if(other!=h&&cfg.hole[other].mass>0)
    dh[i][j]=dh[i][j]+f[other]*F[other]*(s[other].metric[i][j]-Jet(i==j?1:0));
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
  long double Dh[3][3][3]={},DC[3][3][3]={};
  for(int d=0;d<3;d++)for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   Dh[d][i][j]=dh[i][j].d[d+1];
   for(int k=0;k<3;k++)Dh[d][i][j]-=SC[k][d][i].v*dh[k][j].v+SC[k][d][j].v*dh[i][k].v;
  }
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)
   for(int l=0;l<3;l++)DC[i][j][k]+=.5L*b.inv[i][l].v*(Dh[j][k][l]+Dh[k][j][l]-Dh[l][j][k]);
  for(int i=0;i<3;i++){
   for(int j=0;j<3;j++){
    div[i]+=2.0L/3*std::pow(s[h].psi.v,6)*si[i][j].v*s[h].K.d[j+1];
    div[i]+=du[i][j].d[j+1];
    for(int k=0;k<3;k++)div[i]+=SC[i][j][k].v*du[k][j].v+SC[j][j][k].v*du[i][k].v
      +DC[i][j][k]*up[k][j].v+DC[j][j][k]*up[i][k].v;
   }
  }
 }
 for(int i=0;i<3;i++){
  for(int j=0;j<3;j++)div[i]-=b.inv[i][j].v*trace.d[j+1]/3;
  out[i]=(double)div[i];
 }
}

void background(const HiSpID_Config&cfg,const double*x,Background&b){
 Seed s[2];Jet radius[2],F[2]={1,1},f[2]={1,1};
 b.psi=1;b.g=1;b.K=0;
 for(int h=0;h<2;h++)if(cfg.hole[h].mass>0){
  Jet r2=0;for(int i=0;i<3;i++){
   Jet dx=Jet::variable(x[i]-cfg.hole[h].center[i],i+1);r2=r2+dx*dx;
  }radius[h]=sqrt(r2);
  seed(cfg.hole[h],cfg.conformal_choice,x,s[h]);
  if(cfg.far_radius>0)F[h]=exp(-power(radius[h]/Jet(cfg.far_radius),4));
  b.psi=b.psi+F[h]*(s[h].psi-Jet(1));
  b.far_correction=b.far_correction+(Jet(1)-F[h])*(s[h].psi-Jet(1));
  b.g=b.g*inner(radius[h],cfg.inner_min[h],cfg.inner_max[h]);
 }
 for(int h=0;h<2;h++)if(cfg.omega[h]>0&&cfg.hole[1-h].mass>0)
  f[h]=Jet(1)-exp(-power(radius[1-h]/Jet(cfg.omega[h]),cfg.attenuation_power));
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  b.metric[i][j]=Jet(i==j?1:0);b.M[i][j]=0;
  for(int h=0;h<2;h++)if(cfg.hole[h].mass>0){
   b.metric[i][j]=b.metric[i][j]+f[h]*F[h]*(s[h].metric[i][j]-Jet(i==j?1:0));
   b.M[i][j]=b.M[i][j]+s[h].A[i][j];
  }
 }
 for(int h=0;h<2;h++)if(cfg.hole[h].mass>0)b.K=b.K+f[h]*F[h]*s[h].K;
 invert(b.metric,b.inv);connection(b.metric,b.inv,b.C);
 Jet trace;seed_sum_source(cfg,s,f,F,b,trace,b.divM);
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)b.M[i][j]=b.M[i][j]-b.metric[i][j]*trace/Jet(3);
 b.R=(double)curvature(b.inv,b.C);
 b.lapPsi=laplacian(b.inv,b.C,b.psi);
 if(cfg.inner_flatten){
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)
   b.opmetric[i][j]=Jet(i==j?1:0)+b.g*(b.metric[i][j]-Jet(i==j?1:0));
  invert(b.opmetric,b.opinv);
  for(int k=0;k<3;k++)for(int i=0;i<3;i++)for(int j=0;j<3;j++)b.opC[k][i][j]=b.g*b.C[k][i][j];
 }else{
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){b.opmetric[i][j]=b.metric[i][j];b.opinv[i][j]=b.inv[i][j];
   for(int k=0;k<3;k++)b.opC[k][i][j]=b.C[k][i][j];}
 }
 if(!std::isfinite(b.R)||!std::isfinite(b.lapPsi))throw std::runtime_error("nonfinite background geometry");
}

double laplacian(const Jet inv[3][3],const Jet C[3][3][3],const Jet&u){
 long double out=0;for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  long double h=u.h[i+1][j+1];for(int k=0;k<3;k++)h-=C[k][i][j].v*u.d[k+1];
  out+=inv[i][j].v*h;
 }
 return (double)out;
}

void longitudinal(const Jet metric[3][3],const Jet C[3][3][3],const Jet b[3],Jet L[3][3]){
 Jet db[3][3],div=0;for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  db[i][j]=diff(b[j],i+1);for(int k=0;k<3;k++)db[i][j]=db[i][j]+C[j][i][k]*b[k];
  if(i==j)div=div+db[i][j];
 }
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  L[i][j]=-Jet(2.0L/3)*metric[i][j]*div;
  for(int k=0;k<3;k++)L[i][j]=L[i][j]+metric[i][k]*db[j][k]+metric[j][k]*db[i][k];
 }
}

void divergence(const Jet inv[3][3],const Jet C[3][3][3],const Jet T[3][3],double out[3]){
 Jet up[3][3];for(int i=0;i<3;i++)for(int j=0;j<3;j++)
  for(int k=0;k<3;k++)for(int l=0;l<3;l++)up[i][j]=up[i][j]+inv[i][k]*inv[j][l]*T[k][l];
 for(int i=0;i<3;i++){
  long double v=0;for(int j=0;j<3;j++){
   v+=up[i][j].d[j+1];for(int k=0;k<3;k++)v+=C[i][j][k].v*up[k][j].v+C[j][j][k].v*up[i][k].v;
  }out[i]=(double)v;
 }
}

static long double contract(const Background&b,const Jet A[3][3],const Jet D[3][3]){
 long double q=0;for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)for(int l=0;l<3;l++)
  q+=b.inv[i][k].v*b.inv[j][l].v*A[i][j].v*D[k][l].v;
 return q;
}

void equations(const Background&b,const Jet u[4],double out[4],const Jet *du){
 Jet L[3][3],A[3][3];longitudinal(b.metric,b.C,u+1,L);
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)A[i][j]=b.M[i][j]+L[i][j];
 long double psi=b.psi.v+u[0].v;
 if(!(psi>0)){for(int k=0;k<4;k++)out[k]=std::numeric_limits<double>::quiet_NaN();return;}
 long double A2=contract(b,A,A),K=b.K.v,g=b.g.v;
 if(!du){
  out[0]=laplacian(b.opinv,b.opC,u[0])+g*(-psi*b.R/8-std::pow(psi,5)*K*K/12+A2/(8*std::pow(psi,7))+b.lapPsi);
  Jet OL[3][3];longitudinal(b.opmetric,b.opC,u+1,OL);
  divergence(b.opinv,b.opC,OL,out+1);
  for(int i=0;i<3;i++){
   long double grad=0;for(int j=0;j<3;j++)grad+=b.inv[i][j].v*b.K.d[j+1];
   out[i+1]+=g*(b.divM[i]-2.0L/3*std::pow(psi,6)*grad);
  }
 }else{
  Jet DL[3][3],OL[3][3];longitudinal(b.metric,b.C,du+1,DL);
  out[0]=laplacian(b.opinv,b.opC,du[0])+g*((-b.R/8-5*std::pow(psi,4)*K*K/12-7*A2/(8*std::pow(psi,8)))*du[0].v
     +contract(b,A,DL)/(4*std::pow(psi,7)));
  longitudinal(b.opmetric,b.opC,du+1,OL);divergence(b.opinv,b.opC,OL,out+1);
  for(int i=0;i<3;i++){
   long double grad=0;for(int j=0;j<3;j++)grad+=b.inv[i][j].v*b.K.d[j+1];
   out[i+1]-=g*4*std::pow(psi,5)*du[0].v*grad;
  }
 }
}

bool valid(const HiSpID_Config&c,bool sampler_only){
 if(c.conformal_choice<0||c.conformal_choice>1||c.attenuation_power<2||c.attenuation_power%2)return false;
 if(!std::isfinite(c.far_radius)||!std::isfinite(c.tolerance)||c.tolerance<=0||c.max_newton<0||c.max_krylov<1||c.krylov_restart<2||c.krylov_restart>200)return false;
 for(int k=0;k<3;k++)if(c.n[k]<4||c.n[k]>256||(k==2&&c.n[k]%2))return false;
 /* Bound the compact cache, four-field modal FD stencil and Krylov basis before
  * allocation. Larger grids require an explicit per-context budget. */
 if(c.memory_limit_mib<16||c.memory_limit_mib>8192)return false;
 const double npt=(double)c.n[0]*c.n[1]*c.n[2];
 // Sampler: values, coefficients and SpecCoef's two padded scalar scratch
 // arrays/pointer tables, with headroom for the smallest allowed grids.
 const double bytes_per_point=sampler_only?128:3192+4*(5*12+48)+32*(2*c.krylov_restart+30);
 // Two dense matrices per radial block, scalar/vector factor groups, shared
 // by real-Fourier partners. Add permutations/diagonal couplings explicitly.
 const double block_bytes=sampler_only?0:2.0*(c.n[2]/2+1)*c.n[1]*(16.0*c.n[0]*c.n[0]+24.0*c.n[0]);
 if(npt*bytes_per_point+block_bytes>(double)c.memory_limit_mib*1024*1024)return false;
 bool active=false;for(int h=0;h<2;h++){
  const auto&v=c.hole[h];if(!std::isfinite(v.mass)||v.mass<0)return false;
  for(int k=0;k<3;k++)if(!std::isfinite(v.center[k])||!std::isfinite(v.spin[k])||!std::isfinite(v.velocity[k]))return false;
  if(v.mass>0){active=true;if(norm(v.spin)>=v.mass*v.mass||norm(v.velocity)>=1)return false;}
  if(!std::isfinite(c.omega[h])||!std::isfinite(c.inner_min[h])||!std::isfinite(c.inner_max[h]))return false;
  if(c.inner_max[h]>0&&(c.inner_min[h]<0||c.inner_min[h]>=c.inner_max[h]))return false;
 }return active;
}
}

extern "C" {
const char *HiSpID_last_error(void){return hispid::last_error.c_str();}
void HiSpID_default_config(HiSpID_Config*c){
 if(!c)return;std::memset(c,0,sizeof(*c));
 c->hole[0].mass=c->hole[1].mass=.5;c->hole[0].center[0]=3;c->hole[1].center[0]=-3;
 c->n[0]=c->n[1]=16;c->n[2]=8;c->conformal_choice=1;c->inner_flatten=1;
 c->omega[0]=c->omega[1]=.5;c->attenuation_power=4;
 for(int h=0;h<2;h++){c->inner_min[h]=.05;c->inner_max[h]=.1;}
 c->far_radius=40;c->tolerance=1e-10;c->max_newton=12;c->max_krylov=600;c->krylov_restart=40;
 c->memory_limit_mib=2048;
}
int HiSpID_seed(const HiSpID_Hole*h,int choice,int count,const double*xyz,HiSpID_Point*out){
 if(!h||!xyz||!out||count<0||!(h->mass>0)||choice<0||choice>1)return -1;
 HiSpID_Config c;HiSpID_default_config(&c);c.hole[0]=*h;c.hole[1].mass=0;
 if(!hispid::valid(c))return -1;
 try{for(int p=0;p<count;p++){
  hispid::Seed s;hispid::seed(*h,choice,xyz+3*p,s);std::memset(out+p,0,sizeof(*out));
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   out[p].gamma[3*i+j]=s.physical[i][j].v;out[p].Kij[3*i+j]=s.extrinsic[i][j].v;
   out[p].conformal_metric[3*i+j]=s.metric[i][j].v;out[p].Atilde[3*i+j]=s.A[i][j].v;
   }
   out[p].psi=s.psi.v;out[p].mean_curvature=s.K.v;out[p].attenuation=1;
 }}catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
int HiSpID_operators(const HiSpID_Config*c,const double*x,const double*j,double*out){
 if(!c||!x||!j||!out||!hispid::valid(*c))return -1;
 try{
  hispid::Background b;hispid::background(*c,x,b);hispid::Jet u[4];
  const int ii[6]={1,1,1,2,2,3},jj[6]={1,2,3,2,3,3};
  for(int k=0;k<4;k++){u[k].v=j[10*k];for(int d=0;d<3;d++)u[k].d[d+1]=j[10*k+1+d];
   for(int d=0;d<6;d++)u[k].h[ii[d]][jj[d]]=u[k].h[jj[d]][ii[d]]=j[10*k+4+d];}
  out[0]=hispid::laplacian(b.inv,b.C,u[0]);hispid::Jet L[3][3];
  hispid::longitudinal(b.metric,b.C,u+1,L);hispid::divergence(b.inv,b.C,L,out+1);out[4]=b.R;
 }catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
}
