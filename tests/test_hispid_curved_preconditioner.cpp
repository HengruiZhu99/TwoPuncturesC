#include "../src/HiSpID_curved_preconditioner.hpp"
#include "../src/HiSpID_modal_block.hpp"
#include <cstdio>
int main(){
 hispid::Cached flat{};
 for(int d:{4,7,9})flat.lap[d]=1;
 for(int v=0;v<3;v++)for(int d:{4,7,9})flat.vec[v][v][d]=1+(hispid::hessian_i(d-4)==v+1?1./3:0);
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)if(i!=j)for(int d=4;d<10;d++)
  if((hispid::hessian_i(d-4)==i+1&&hispid::hessian_j(d-4)==j+1)||(hispid::hessian_i(d-4)==j+1&&hispid::hessian_j(d-4)==i+1))flat.vec[i][j][d]=1./3;
 double worst=0;
 for(double sigma:{.0000018592359,.0002,.008,.13,.65,.97,.9999981407641})for(double zeta:{-.9999247,-.995,-.71,.13,.97,.9999247})for(int m:{0,1,2,3,4,7,8})for(int group:{0,1}){
  double b=8,t=hispid::AxisDerivatives::radial_t(sigma),ts,tss;hispid::AxisDerivatives::radial_forward(sigma,ts,tss);
  double k=hispid::AxisDerivatives::angular_stretch,T=std::tanh(k),e=std::tanh(k*zeta)/T,s2=1-e*e,C=(1+t)/(1-t),h2=C*C-1,D=b*b*(h2+s2),rho=b*std::sqrt(h2*s2),phi=.37;
  double de=k*(1-T*T*e*e)/T,dde=-2*k*k*e*(1-T*T*e*e),r=hispid::AxisDerivatives::exponent(m),mu=group&&m!=1?10./9:1;
  auto q=hispid::curved_symbol(flat,b,b*C*e,rho*std::cos(phi),rho*std::sin(phi),m,group);
  double ctt=t*(1-t)*(1-t),ct=(1-t)*((r+1)*(1-t)-2*t);
  double expected[5]={mu*(-(r+1)*(r+1-t)+(r*r-m*m)*(1/h2+1/s2))/D,
    mu*(ct/ts-ctt*tss/(ts*ts*ts))/D,mu*(-2*(r+1)*e/de-s2*dde/(de*de*de))/D,
    mu*ctt/(ts*ts*D),mu*s2/(de*de*D)};
  worst=std::max(worst,std::abs(q.mixed)/(1+2*std::sqrt(q.second[0]*q.second[1])));
  double got[5]={q.zero,q.first[0],q.first[1],q.second[0],q.second[1]};
  for(int d=0;d<5;d++){
   if(!std::isfinite(got[d]))return 1;
   worst=std::max(worst,std::abs(got[d]-expected[d])/(1+std::abs(expected[d])));
  }
 }
 printf("curved modal symbol independent flat scalar/vector limit: %.17g\n",worst);
 if(worst>=1e-9)return 1;
 // Axial rotation is a Killing field of the stationary conformal Kerr
 // metric. Its m=1 pair has P proportional to D^2, D=b(C+1). This tests
 // the connection/component cancellation, independently of any inverse.
 HiSpID_Config axis;HiSpID_default_config(&axis);const double focus=.5;
 axis.hole[0]={};axis.hole[0].mass=1;axis.hole[0].center[0]=focus;axis.hole[0].spin[0]=.9999;
 axis.hole[1]={};axis.hole[1].center[0]=-focus;axis.far_radius=0;axis.conformal_choice=1;
 axis.inner_min[0]=.1;axis.inner_max[0]=.2;
 double killing_error=0,mean_killing_defect=0;
 for(auto ray:std::array<std::array<double,2>,4>{{{{0,.3}},{{.35,.12}},{{.8,.2}},{{-.3,.25}}}}){
  const double phi=.37,x[3]={focus+ray[0],ray[1]*std::cos(phi),ray[1]*std::sin(phi)};
  hispid::BackgroundT<long double>bg;if(hispid::background_geometry(axis,x,bg,true,HISPID_SEED_LIU))return 1;
  hispid::Cached c;hispid::cache(bg,c);
  double plus=std::hypot(x[0]+focus,ray[1]),minus=std::hypot(x[0]-focus,ray[1]);
  double C=(plus+minus)/(2*focus),t=(C-1)/(C+1),sigma=hispid::AxisDerivatives::radial_inverse(t),ts,tss;
  hispid::AxisDerivatives::radial_forward(sigma,ts,tss);
  auto q=hispid::curved_symbol(c,focus,x[0],x[1],x[2],1,1);
  double d=2*ts/(1-t),dd=2*tss/(1-t)+6*ts*ts/((1-t)*(1-t));
  double terms[3]={q.zero,q.first[0]*d,q.second[0]*dd};
  killing_error=std::max(killing_error,std::abs(terms[0]+terms[1]+terms[2])/(1+std::abs(terms[0])+std::abs(terms[1])+std::abs(terms[2])));
  auto mean=hispid::curved_symbol(c,focus,x[0],x[1],x[2],1,0);
  double old_terms[3]={mean.zero,mean.first[0]*d,mean.second[0]*dd};
  mean_killing_defect=std::max(mean_killing_defect,std::abs(old_terms[0]+old_terms[1]+old_terms[2])/(1+std::abs(old_terms[0])+std::abs(old_terms[1])+std::abs(old_terms[2])));
 }
 printf("stationary Kerr axial Killing symbol witness: %.17g\n",killing_error);
 printf("scalar-only axial Killing symbol defect: %.17g\n",mean_killing_defect);
 if(killing_error>=1e-9)return 1;
 // Independent full matrix witness for the three-band polar couplings.
 const int na=7,nb=11,N=na*nb;
 hispid::ModalBlock B;B.na=na;B.nb=nb;B.banded_coupling=true;
 B.lu.assign(na*na*nb,0);B.transfer.resize(na*na*nb);B.lower.assign(3*N,0);B.upper.assign(3*N,0);B.permutation.resize(N);
 std::vector<double>A(N*N,0),exact(N),rhs(4*N),multiple(2*N);
 for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
  int row=j*na+i;exact[row]=std::sin(.37*row);
  for(int k=0;k<na;k++)A[row*N+j*na+k]=B.lu[row*na+k]=(i==k?6.:.03*std::cos(i+3*k));
  for(int k=std::max(0,i-1);k<=std::min(na-1,i+1);k++){
   if(j)A[row*N+(j-1)*na+k]=B.lower[3*row+k-i+1]=.2+.05*(i-k);
   if(j+1<nb)A[row*N+(j+1)*na+k]=B.upper[3*row+k-i+1]=-.1+.04*(i-k);
  }
 }
 for(int i=0;i<N;i++){double value=0;for(int j=0;j<N;j++)value+=A[i*N+j]*exact[j];rhs[4*i]=value;multiple[2*i]=value;multiple[2*i+1]=2*value;}
 B.factor();B.solve(rhs.data(),0,0);B.solve_many(multiple.data(),2);
 double inverse_error=0;
 for(int i=0;i<N;i++){inverse_error=std::max(inverse_error,std::abs(rhs[4*i]-exact[i]));for(int k=0;k<2;k++)inverse_error=std::max(inverse_error,std::abs(multiple[2*i+k]-(k+1)*exact[i]));}
 printf("mixed-derivative banded inverse full-matrix witness: %.17g\n",inverse_error);
 if(inverse_error>=1e-12)return 1;
 // The requested radially resolved sequence exceeds the previous 512-node
 // cap. Exercise the larger extent/budget with a small sampling-only grid.
 HiSpID_Config c;HiSpID_default_config(&c);c.n[0]=576;c.n[1]=c.n[2]=4;c.memory_limit_mib=HISPID_MAX_MEMORY_MIB;
 auto*sampler=HiSpID_create_sampler(&c);if(!sampler)return 1;HiSpID_destroy(sampler);
 c.n[0]=HISPID_MAX_RADIAL_POINTS+1;if(HiSpID_create_sampler(&c))return 1;
 c.n[0]=4;c.memory_limit_mib=HISPID_MAX_MEMORY_MIB+1;if(HiSpID_create_sampler(&c))return 1;
 return 0;
}
