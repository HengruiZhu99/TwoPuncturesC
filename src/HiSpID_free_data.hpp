#ifndef HISPID_FREE_DATA_HPP
#define HISPID_FREE_DATA_HPP
#include "HiSpID_internal.hpp"
namespace hispid {
inline bool valid_wave22(const HiSpID_Wave22&w){
 const double a[]={w.amplitude,w.omega,w.cycles,w.r0,w.time,w.cutoff_min,w.cutoff_max};
 for(double x:a)if(!std::isfinite(x))return false;
 return w.struct_size==sizeof(w)&&(w.propagation==1||w.propagation==-1)&&
  (w.handedness==1||w.handedness==-1)&&w.omega!=0&&w.cycles>0&&w.cutoff_min>0&&
  w.cutoff_max>w.cutoff_min&&std::isfinite(2*std::acos(-1.)*w.cycles/std::abs(w.omega));
}
/* Pole-regular Cartesian TT projection of a rotating constant polarization.
 * Equivalent to SpEC's two spin-weighted harmonics, without atan2 or 1/rho. */
inline void wave22(const HiSpID_Wave22&w,const double*center,const double*point,
                   Jet h[3][3],Jet ht[3][3]){
 if(w.amplitude==0)return;
 Jet x[3],r2;for(int i=0;i<3;i++){x[i]=Jet::variable(point[i]-center[i],i+1);r2=r2+x[i]*x[i];}
 if(r2.v<=w.cutoff_min*w.cutoff_min)return;
 Jet r=sqrt(r2),cutoff=inner(r,w.cutoff_min,w.cutoff_max);
 Jet y=r-Jet(w.r0)+Jet(w.propagation)*Jet::variable(w.time,0);
 Jet end=Jet(2*std::acos(-1.L)*w.cycles/std::abs(w.omega))-y;
 Jet e0=exp(Jet(-2)*(y.v<0?-y:y)),e1=exp(Jet(-2)*(end.v<0?-end:end));
 Jet f0=e0/(Jet(1)+e0),f1=e1/(Jet(1)+e1);
 Jet window=y.v<0?Jet(2)*(f0-f1):(end.v<0?Jet(2)*(f1-f0):Jet(2)*(Jet(1)-f0-f1));
 Jet angle=Jet(w.handedness)*(Jet(w.omega)*y+Jet(std::acos(-1.L)/4));
 Jet cc=cos(angle),ss=sin(angle),P[3][3],E[3][3],PEP[3][3],trace;
 E[0][0]=cc;E[1][1]=-cc;E[0][1]=E[1][0]=ss;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)P[i][j]=Jet(i==j?1:0)-x[i]*x[j]/r2;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  trace=trace+P[i][j]*E[i][j];
  for(int k=0;k<3;k++)for(int l=0;l<3;l++)PEP[i][j]=PEP[i][j]+P[i][k]*E[k][l]*P[l][j];
 }
 Jet scale=Jet(std::sqrt(15/std::acos(-1.L))*w.amplitude)*cutoff*window/r;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  h[i][j]=scale*(PEP[i][j]-Jet(.5L)*P[i][j]*trace);ht[i][j]=diff(h[i][j],0);
 }
}
/* Only called outside the wave-free throat neighbourhood. No change to the
 * native cancellation-free Kerr K at, or on either side of, the throat. */
inline void stationary_gauge(const HiSpID_Hole&h,const double*point,int family,Jet&alpha,Jet beta[3]){
 const long double m=h.mass,smag=point_norm<long double>(h.spin),a=smag/m;
 long double axis[3]={0,0,1};if(smag>0)for(int i=0;i<3;i++)axis[i]=h.spin[i]/smag;
 Jet x[3],r2,z;for(int i=0;i<3;i++){x[i]=Jet::variable(point[i]-h.center[i],i+1);r2=r2+x[i]*x[i];z=z+Jet(axis[i])*x[i];}
 Jet r=sqrt(r2),mu=z/r;const long double rp=m+std::sqrt(m*m-a*a),rm=a*a/rp,c=rp/4;
 const bool liu=family==HISPID_SEED_LIU;
 Jet rb=liu?r+Jet(2*c)+Jet(c*c)/r:r+Jet(m)+Jet((m*m-a*a)/4)/r;
 Jet radial=liu?(r-Jet(c))*sqrt((rb-Jet(rm))/r):r-Jet((m*m-a*a)/4)/r;
 Jet sigma=rb*rb+Jet(a*a)*mu*mu,AA=power(rb*rb+Jet(a*a),2)-Jet(a*a)*radial*radial*(Jet(1)-mu*mu);
 alpha=radial*sqrt(sigma/AA);
 for(int i=0;i<3;i++)beta[i]=-Jet(2*m*a)*rb/AA*(Jet(axis[(i+1)%3])*x[(i+2)%3]-Jet(axis[(i+2)%3])*x[(i+1)%3]);
}
/* delta_gamma/dt_delta_gamma are physical lab-frame perturbations pulled
 * back to the solver frame. psi_seed and its time independence are fixed. */
inline void perturb_background(const HiSpID_Config&cfg,int family,const double*x,
                                const Jet h[3][3],const Jet ht[3][3],Background&b){
 bool nonzero=false;for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  nonzero|=h[i][j].v!=0||ht[i][j].v!=0;
  for(int d=1;d<4;d++){nonzero|=h[i][j].d[d]!=0||ht[i][j].d[d]!=0;for(int e=1;e<4;e++)nonzero|=h[i][j].h[d][e]!=0;}
 }
 if(!nonzero)return;
 const HiSpID_Hole&hole=cfg.hole[cfg.hole[0].mass>0?0:1];
 Seed s;int status=seed_geometry(hole,cfg.conformal_choice,x,s,true,family);
 if(status)throw std::runtime_error("invalid free-data reference seed");
 Jet alpha,beta[3];stationary_gauge(hole,x,family,alpha,beta);
 if(!(alpha.v>1e-6L))throw std::runtime_error("free-data perturbation must vanish inside and near the Kerr throat");
 Jet psi4=power(s.psi,4),Kij[3][3];b.K=0;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  b.metric[i][j]=s.metric[i][j]+h[i][j]/psi4;
  Jet dt=ht[i][j];for(int k=0;k<3;k++)dt=dt-beta[k]*diff(h[i][j],k+1)-h[k][j]*diff(beta[k],i+1)-h[i][k]*diff(beta[k],j+1);
  Kij[i][j]=s.extrinsic[i][j]-dt/(Jet(2)*alpha);
 }
 if(!(b.metric[0][0].v>0)||!((b.metric[0][0]*b.metric[1][1]-b.metric[0][1]*b.metric[0][1]).v>0)||!invert_checked(b.metric,b.inv))
  throw std::runtime_error("perturbed metric is not positive definite");
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)b.K=b.K+b.inv[i][j]*Kij[i][j]/psi4;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)b.M[i][j]=s.psi*s.psi*(Kij[i][j]-psi4*b.metric[i][j]*b.K/Jet(3));
 connection(b.metric,b.inv,b.C);b.R=curvature(b.inv,b.C);b.lapPsi=laplacian(b.inv,b.C,b.psi);
 // Exact seed identities apply only to the unperturbed seed. Recompute all
 // new source terms from the perturbed tensors, including their derivatives.
 divergence(b.inv,b.C,b.M,b.divM);
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  b.opmetric[i][j]=cfg.inner_flatten?Jet(i==j?1:0)+b.g*(b.metric[i][j]-Jet(i==j?1:0)):b.metric[i][j];
  for(int k=0;k<3;k++)b.opC[k][i][j]=cfg.inner_flatten?b.g*b.C[k][i][j]:b.C[k][i][j];
 }
 invert(b.opmetric,b.opinv);
 if(!std::isfinite(b.R)||!std::isfinite(b.lapPsi)||!std::isfinite(b.K.v))throw std::runtime_error("nonfinite perturbed geometry");
 for(double d:b.divM)if(!std::isfinite(d))throw std::runtime_error("nonfinite perturbed momentum source");
}
}
#endif
