// Test-only access to private operators; no extra public API entry points.
#include "../src/HiSpID_solver.cpp"
#include <iostream>
#include <complex>
// Independent continuous Cartesian expression for the complete retained
// polynomial, using puncture distances and all six Hessians. No legacy
// A/B transform or differentiation matrix is reused by this oracle.
Fields cartesian_polynomial(const HiSpID_Data&s,const double*x){
 using hispid::Jet;Jet X=Jet::variable(x[0],1),Y=Jet::variable(x[1],2),Z=Jet::variable(x[2],3);
 Jet rp=sqrt((X-Jet(s.b))*(X-Jet(s.b))+Y*Y+Z*Z),rm=sqrt((X+Jet(s.b))*(X+Jet(s.b))+Y*Y+Z*Z);
 Jet D=(rp+rm+Jet(2*s.b))/Jet(2),t=Jet(1)-Jet(2*s.b)/D,eta=(rm-rp)/Jet(2*s.b);
 const double lambda=hispid::AxisDerivatives::radial_stretch,kappa=hispid::AxisDerivatives::angular_stretch;
 Jet z=Jet(2)*t/(Jet(lambda)+Jet(1-lambda)*t)-Jet(1),angle=eta*Jet(std::tanh(kappa));
 Jet zeta=hispid::unary(angle,std::atanh(angle.v),1/(1-angle.v*angle.v),2*angle.v/std::pow(1-angle.v*angle.v,2))/Jet(kappa);
 const int na=s.local.n[0],nb=s.local.n[1],np=s.local.n[2],half=np/2;
 std::vector<Jet>Ta(na),Tb(nb);Ta[0]=Tb[0]=1;Ta[1]=z;Tb[1]=zeta;
 for(int i=2;i<na;i++)Ta[i]=Jet(2)*z*Ta[i-1]-Ta[i-2];
 for(int j=2;j<nb;j++)Tb[j]=Jet(2)*zeta*Tb[j-1]-Tb[j-2];
 Ta[0]=Tb[0]=.5;Jet result[4];
 for(int k=0;k<np;k++){
  int m=k<=half?k:k-half,r=hispid::AxisDerivatives::exponent(m);Jet re(1),im(0);
  for(int q=0;q<m;q++){Jet next=re*Y-im*Z;im=re*Z+im*Y;re=next;}
  Jet G=k<=half?re:im;if(r!=m)G=G*power(Y*Y+Z*Z,.5L*(r-m));
  Jet factor=Jet(-4*s.b*std::sqrt((m==0||m==half?1.:2.)/np))*G/power(D,r+1);
  for(int v=0;v<4;v++){
   Jet P;for(int j=0;j<nb;j++)for(int i=0;i<na;i++)P=P+Jet(s.coefficients[v+4*(i+na*(j+nb*k))])*Ta[i]*Tb[j];
   result[v]=result[v]+factor*P;
  }
 }
 Fields out{};for(int v=0;v<4;v++){
  out[v][0]=result[v].v;for(int d=1;d<4;d++)out[v][d]=result[v].d[d];
  for(int d=0;d<6;d++)out[v][4+d]=result[v].h[hi[d]][hj[d]];
 }return out;
}
// A coupled Navier inverse control complements the scalar Laplace controls.
// RHS Hessians come directly from the Cartesian distance expression.
bool flat_vector_inverse(bool anisotropic=false,int restart=32,int max_krylov=300,bool difficult_only=false){
 bool passed=true;
 for(int N:difficult_only?std::vector<int>{32}:std::vector<int>{8,16,32}){
  HiSpID_Config cfg;HiSpID_default_config(&cfg);cfg.n[0]=cfg.n[1]=N;cfg.n[2]=16;
  cfg.conformal_choice=0;cfg.inner_max[0]=cfg.inner_max[1]=0;cfg.far_radius=0;
  HiSpID_Data*s=HiSpID_create(&cfg);if(!s)return false;s->basefields.resize(s->npt);
  for(double epsilon:difficult_only?std::vector<double>{.4}:anisotropic?std::vector<double>{.01,.4}:std::vector<double>{0.}){
  hispid::Background flat{};flat.psi=1;double metric_inverse[3][3];const double direction[3]={.5,.8,.3},length2=.25+.64+.09;
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   flat.metric[i][j]=flat.opmetric[i][j]=(i==j?1.:0.)+epsilon*direction[i]*direction[j]/length2;
   metric_inverse[i][j]=(i==j?1.:0.)-epsilon/(1+epsilon)*direction[i]*direction[j]/length2;
   flat.inv[i][j]=flat.opinv[i][j]=metric_inverse[i][j];
  }
  for(auto&g:s->geometry){double weight=g.weight;cache(flat,g);g.weight=weight;}
  Sparse M=preconditioner(*s);s->local.krylov_restart=restart;s->local.max_krylov=max_krylov;
  for(int mode:difficult_only?std::vector<int>{0}:std::vector<int>{0,1,2,4})for(int component:difficult_only?std::vector<int>{1,2}:anisotropic?std::vector<int>{0,1,2}:std::vector<int>{1,2}){
   const double amplitude=1e-4,normal=std::sqrt((mode?2.:1.)/16);
   auto reference=[&](const double*x){
    using hispid::Jet;Jet X=Jet::variable(x[0],1),Y=Jet::variable(x[1],2),Z=Jet::variable(x[2],3);
    Jet rp=sqrt((X-Jet(s->b))*(X-Jet(s->b))+Y*Y+Z*Z),rm=sqrt((X+Jet(s->b))*(X+Jet(s->b))+Y*Y+Z*Z),D=(rp+rm+Jet(2*s->b))/Jet(2);
    Jet re(1),im(0);for(int q=0;q<mode;q++){Jet next=re*Y-im*Z;im=re*Z+im*Y;re=next;}
    return Jet(-4*s->b*amplitude*normal)*re/power(D,mode+1);
   };
   std::vector<double>rhs(s->ntotal,0),known(s->ntotal,0),action(s->ntotal),difference(s->ntotal),solution;
   for(int k=0;k<16;k++)for(int j=0;j<N;j++)for(int i=0;i<N;i++){
    int p=pindex(*s,i,j,k);double x[3];Fields zero{};transform(*s,i,j,k,zero,x);auto v=reference(x);
    double trace=0;for(int d=0;d<3;d++)for(int e=0;e<3;e++)trace+=metric_inverse[d][e]*v.h[d+1][e+1];
    if(component==0)rhs[4*p]=s->geometry[p].weight*trace;
    else for(int d=0;d<3;d++){
     double value=d+1==component?trace:0;
     for(int e=0;e<3;e++)value+=metric_inverse[d][e]*v.h[e+1][component]/3;
     rhs[4*p+d+1]=s->geometry[p].weight*value;
    }
   }
   for(int p=mode*N*N;p<(mode+1)*N*N;p++)known[4*p+component]=amplitude;
   jvp(*s,known.data(),action.data());for(int p=0;p<s->ntotal;p++)difference[p]=rhs[p]-action[p];
   const double exact_floor=norm2v(difference)/norm2v(rhs);const int before=s->diag.krylov_iterations;
   bool converged=gmres(*s,M,rhs,solution,1e-11);
   jvp(*s,solution.data(),action.data());for(int p=0;p<s->ntotal;p++)difference[p]=rhs[p]-action[p];
   const double final_error=norm2v(difference)/norm2v(rhs);s->values=solution;s->coefficients_valid=false;
   double field_error=0;const double points[3][3]={{.7,.8,.4},{s->b+.003,.005,.002},{-s->b+.003,.006,.001}};
   for(const auto&x:points){Fields recovered{};sample_fields(*s,x,recovered);auto v=reference(x);
    for(int c=0;c<4;c++)for(int d=0;d<4;d++){
     const double expected=c==component?(d?double(v.d[d]):double(v.v)):0;
     field_error=std::max(field_error,std::abs(recovered[c][d]-expected));
    }
   }
   if(anisotropic)std::cout<<"constant anisotropic metric epsilon"<<epsilon<<" ";
   std::cout<<"exact delta coupled vector inverse N"<<N<<" mode"<<mode<<" component"<<component<<": converged "<<converged<<", Krylov "<<s->diag.krylov_iterations-before<<", exact-P linear floor "<<exact_floor<<", final linear relative residual "<<final_error<<", field/gradient error "<<field_error<<'\n';
   passed=passed&&converged&&field_error<1e-11;
  }
  }HiSpID_destroy(s);
 }return passed;
}
int main(int argc,char**argv){
 std::cout<<std::unitbuf;
 if(argc>1&&std::string(argv[1])=="--vector-only")return flat_vector_inverse()?0:1;
 if(argc>1&&std::string(argv[1])=="--anisotropic-vector-only")return flat_vector_inverse(true)?0:1;
 if(argc>1&&std::string(argv[1])=="--anisotropic-difficult-only")return flat_vector_inverse(true,128,1024,true)?0:1;
 double gradient_error=0;HiSpID_Config cfg;HiSpID_default_config(&cfg);cfg.n[0]=12;cfg.n[1]=13;cfg.n[2]=16;
 cfg.conformal_choice=0;cfg.inner_max[0]=cfg.inner_max[1]=0;cfg.far_radius=0;
 HiSpID_Data*s=HiSpID_create(&cfg);if(!s)return 1;
 unsigned random=42;for(auto&v:s->values){random=1664525*random+1013904223;v=1e-8*(double(random)/4294967296.-.5);}
 std::vector<Fields>f;fields(*s,s->values.data(),f);make_coefficients(*s);
 for(int k:{0,1,3,7})for(int j:{1,5,10})for(int i:{1,5,10}){
  double x[3];Fields zero{};transform(*s,i,j,k,zero,x);Fields sampled{};sample_fields(*s,x,sampled);
  const int p=pindex(*s,i,j,k);
  for(int v=0;v<4;v++)for(int d=0;d<4;d++)gradient_error=std::max(gradient_error,std::abs(sampled[v][d]-f[p][v][d])/(1+std::abs(f[p][v][d])));
 }
 HiSpID_destroy(s);std::cout<<"random modal fields: collocation/sampler values and Cartesian gradients "<<gradient_error<<'\n';
 bool passed=gradient_error<1e-11;
 for(int N:{12,24,40}){
  cfg.n[0]=N;cfg.n[1]=N+1;s=HiSpID_create(&cfg);if(!s)return 1;
  for(auto&v:s->values){random=1664525*random+1013904223;v=1e-8*(double(random)/4294967296.-.5);}
  fields(*s,s->values.data(),f);make_coefficients(*s);double hessian_error=0;
  for(int k:{0,1,3,7})for(int j:{0,N/2,N})for(int i:{0,N/2,N-1}){
   double point[3];Fields zero{};transform(*s,i,j,k,zero,point);
   const Fields reference=cartesian_polynomial(*s,point);const int p=pindex(*s,i,j,k);
   for(int v=0;v<4;v++)for(int d=0;d<10;d++)hessian_error=std::max(hessian_error,std::abs(reference[v][d]-f[p][v][d])/(1+std::abs(reference[v][d])));
  }
  std::cout<<"random high-degree Cartesian all-Hessian control N"<<N<<": normalized error "<<hessian_error<<'\n';
  passed=passed&&hessian_error<1e-8;HiSpID_destroy(s);
 }
 if(argc>1&&std::string(argv[1])=="--hessian-only")return passed?0:1;
 for(int N:{8,16,32,64}){
  cfg.n[0]=cfg.n[1]=N;cfg.n[2]=16;s=HiSpID_create(&cfg);if(!s)return 1;
  s->basefields.resize(s->npt);for(auto&g:s->geometry){g.g=0;g.R=g.K=0;g.psi=1;std::fill(g.inv,g.inv+9,0);g.inv[0]=g.inv[4]=g.inv[8]=1;}
  Sparse M=preconditioner(*s);std::vector<double>rhs(s->ntotal,0),solution;
  // Independent matrix-times-known-vector check exercises Schur complements,
  // radial/polar boundaries, scalar/vector reuse and sine/cosine reuse.
  std::vector<double>known(s->ntotal),modal_rhs(s->ntotal);
  for(auto&v:known){random=1664525*random+1013904223;v=double(random)/4294967296.-.5;}
  for(int row=0;row<s->ntotal;row++)for(size_t q=0;q<M.col[row].size();q++)modal_rhs[row]+=M.val[row][q]*known[M.col[row][q]];
  for(int k=0;k<cfg.n[2];k++)for(int v=0;v<4;v++)M.blocks[2*(k<=8?k:k-8)+(v?1:0)].solve(modal_rhs,k,v);
  double inverse_error=0;for(int row=0;row<s->ntotal;row++)inverse_error=std::max(inverse_error,std::abs(modal_rhs[row]-known[row]));
  std::cout<<"modal FD direct inverse N"<<N<<": maximum known-vector error "<<inverse_error<<'\n';
  passed=passed&&inverse_error<5e-10;
  for(int k=0;k<cfg.n[2];k++)for(int j=0;j<N;j++)for(int i=0;i<N;i++){
   const int p=pindex(*s,i,j,k);const double a=.5*(s->derivatives.coordinate[0][i]+1),t=a*a,B=s->derivatives.coordinate[1][j],eta=-2*B/(1+B*B),D=s->b*s->b*(4*t/((1-t)*(1-t))+1-eta*eta);
   rhs[4*p]=s->geometry[p].weight*(-2*(1-t))/D*(-1+t)/std::sqrt(double(cfg.n[2]));
  }
  M.solve(rhs,solution);double low=0,high=0;
  for(int p=0;p<s->npt;p++){if(p<N*N)low=std::max(low,std::abs(solution[4*p]-1));else high=std::max(high,std::abs(solution[4*p]));}
  s->values=solution;s->coefficients_valid=false;Fields point{};double x[3]={.7,.8,.4};sample_fields(*s,x,point);
  std::cout<<"preconditioner low-mode leakage N"<<N<<": m0 error "<<low<<", other P max "<<high<<", sampled u "<<point[0][0]<<'\n';
  // The FD preconditioner approximates the spectral operator: m0 need not
  // equal the exact spectral solution even though its FD inverse is direct.
  // Inspect leakage separately without filtering the PDE or accepted data.
  for(int p=0;p<N*N;p++)solution[4*p]=0;
  s->values=solution;s->coefficients_valid=false;Fields leakage{};sample_fields(*s,x,leakage);
  double leak=0;for(int d=0;d<4;d++)leak=std::max(leak,std::abs(leakage[0][d]));
  std::cout<<"  higher-mode reconstructed value/gradient max "<<leak<<'\n';
  passed=passed&&high==0&&leak==0;
  // An analytic flat-Laplacian source, independent of the Cartesian Hessian
  // transform, must recover the prescribed field without spurious modes.
  s->local.krylov_restart=32;s->local.max_krylov=2000;
  const auto rounded_background=s->geometry;
  for(bool exact_flat:{false,true}){
  if(exact_flat){
   // A transform-only control must not inherit quotient/connection rounding
   // from the Brill-Lindquist seed. Build its complete cached operator from
   // exact constant delta jets, retaining only the original row weights.
   hispid::Background flat{};flat.psi=hispid::Jet(1);
   for(int i=0;i<3;i++)flat.metric[i][i]=flat.inv[i][i]=flat.opmetric[i][i]=flat.opinv[i][i]=hispid::Jet(1);
   for(auto&g:s->geometry){const double weight=g.weight;cache(flat,g);g.weight=weight;}
  }
  for(int mode:{0,1,2,4}){
   const double amplitude=1e-4;std::fill(rhs.begin(),rhs.end(),0);
   for(int k=0;k<cfg.n[2];k++)for(int j=0;j<N;j++)for(int i=0;i<N;i++){
    const int p=pindex(*s,i,j,k);const double a=.5*(s->derivatives.coordinate[0][i]+1),t=a*a,B=s->derivatives.coordinate[1][j],eta=-2*B/(1+B*B),sn=std::sqrt(1-eta*eta);
    const double denominator=s->b*s->b*(4*t/((1-t)*(1-t))+1-eta*eta);
    rhs[4*p]=s->geometry[p].weight*(-2*(1-t)*std::pow(a*sn,mode))/denominator
       *(-(mode+1)*(mode+1-t))*amplitude*std::sqrt((mode?2.:1.)/cfg.n[2])*std::cos(mode*2*Pi*k/cfg.n[2]);
   }
   std::vector<double>reference(s->ntotal,0),action(s->ntotal),difference(s->ntotal);
   for(int p=mode*N*N;p<(mode+1)*N*N;p++)reference[4*p]=amplitude;
   jvp(*s,reference.data(),action.data());for(int p=0;p<s->ntotal;p++)difference[p]=rhs[p]-action[p];
   const double exact_floor=norm2v(difference)/norm2v(rhs);const int before=s->diag.krylov_iterations;
   const bool converged=gmres(*s,M,rhs,solution,1e-13);
   jvp(*s,solution.data(),action.data());for(int p=0;p<s->ntotal;p++)difference[p]=rhs[p]-action[p];
   const double final_linear_error=norm2v(difference)/norm2v(rhs);
   s->values=solution;s->coefficients_valid=false;Fields recovered{};sample_fields(*s,x,recovered);
   double dp[3]={x[0]-s->b,x[1],x[2]},dm[3]={x[0]+s->b,x[1],x[2]},rp=0,rm=0;
   for(int d=0;d<3;d++){rp+=dp[d]*dp[d];rm+=dm[d]*dm[d];}rp=std::sqrt(rp);rm=std::sqrt(rm);
   const double D=(rp+rm+2*s->b)/2,pref=-4*s->b*amplitude*std::sqrt((mode?2.:1.)/cfg.n[2]);
   const std::complex<double>w(x[1],x[2]),wm=std::pow(w,mode),Wi[3]={{0,0},{1,0},{0,1}};
   double error=std::abs(recovered[0][0]-pref*std::real(wm)/std::pow(D,mode+1));
   for(int d=0;d<3;d++){
    const double dD=.5*(dp[d]/rp+dm[d]/rm);
    std::complex<double>derivative=-double(mode+1)*wm*dD/std::pow(D,mode+2);
    if(mode)derivative+=double(mode)*std::pow(w,mode-1)*Wi[d]/std::pow(D,mode+1);
    error=std::max(error,std::abs(recovered[0][d+1]-pref*std::real(derivative)));
   }
   double other=0;for(int p=0;p<s->npt;p++)if(p/(N*N)!=mode)other=std::max(other,std::abs(solution[4*p]));
   std::cout<<"analytic "<<(exact_flat?"exact delta":"BL-rounded delta")<<" spectral inverse N"<<N<<" mode"<<mode<<": converged "<<converged<<", Krylov "<<s->diag.krylov_iterations-before<<", exact-P linear floor "<<exact_floor<<", final linear relative residual "<<final_linear_error<<", field/gradient error "<<error<<", other-mode P max "<<other<<'\n';
   passed=passed&&converged&&error<1e-7*amplitude;
  }
  }
  s->geometry=rounded_background;
  if(N==16){
   // Distinguish the scalar potential factor from the shared vector factors
   // using a bounded spatially varying positive curvature coefficient.
   for(int p=0;p<s->npt;p++){
    const int i=p%N,j=(p/N)%N;s->geometry[p].g=1;
    s->geometry[p].R=.2*(1+.1*s->derivatives.coordinate[0][i]+.1*s->derivatives.coordinate[1][j]);
   }
   Sparse distinct=preconditioner(*s);std::fill(modal_rhs.begin(),modal_rhs.end(),0);
   for(int row=0;row<s->ntotal;row++)for(size_t q=0;q<distinct.col[row].size();q++)modal_rhs[row]+=distinct.val[row][q]*known[distinct.col[row][q]];
   for(int k=0;k<cfg.n[2];k++)for(int v=0;v<4;v++)distinct.blocks[2*(k<=8?k:k-8)+(v?1:0)].solve(modal_rhs,k,v);
   double distinct_error=0;for(int row=0;row<s->ntotal;row++)distinct_error=std::max(distinct_error,std::abs(modal_rhs[row]-known[row]));
   std::cout<<"modal scalar potential/vector reuse: maximum known-vector error "<<distinct_error<<'\n';
   passed=passed&&distinct_error<5e-10;
  }
  HiSpID_destroy(s);
 }
 return passed&&flat_vector_inverse()?0:1;
}
