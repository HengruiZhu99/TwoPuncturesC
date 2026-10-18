#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include "../src/HiSpID_cache_kernels.hpp"
using J=hispid::Jet;
static J cheb(int n,J x){J a(1),b=x;if(!n)return a;for(int k=2;k<=n;k++){J c=J(2)*x*b-a;a=b;b=c;}return b;}
int main(){
 double value_error=0,K_error=0;const double pi=std::acos(-1.0),lambda=HISPID_RADIAL_STRETCH,kappa=HISPID_ANGULAR_STRETCH;
 for(int na:{16,192}){
  HiSpID_Config c;HiSpID_default_config(&c);c.n[0]=na;c.n[1]=2*na;c.n[2]=16;c.conformal_choice=0;c.far_radius=0;
  c.hole[0]={.6,{3,0,0},{.072,.054,.108},{.03,.06,.01}};c.hole[1]={.4,{-3,0,0},{-.032,.048,.016},{-.02,-.07,.025}};
  auto*s=HiSpID_create_with_seed_family(&c,HISPID_SEED_TRUMPET_R0_M,0,0,1);if(!s)return 1;
  const int nb=c.n[1],np=c.n[2];std::vector<double>u(size_t(4)*na*nb*np,0.);
  for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
   double P=100+1e-6*std::cos((na-1)*std::acos(-std::cos(pi*(i+.5)/na)))*std::cos((nb-1)*std::acos(-std::cos(pi*(j+.5)/nb)));
   u[4*(i+na*(j+nb*4))+2]=P;
  }
  if(HiSpID_set_unknowns(s,u.data(),u.size()))return 2;
  for(auto z: {std::pair<double,double>{.13,.37},{-.42,-.23}}){
   const double sig=(1+z.first)/2,t=lambda*sig/(1-(1-lambda)*sig),eta=std::tanh(kappa*z.second)/std::tanh(kappa),phi=.31;
   const double rho=6*std::sqrt(t)*std::sqrt(1-eta*eta)/(1-t),x[3]={3*(1+t)*eta/(1-t),rho*std::cos(phi),rho*std::sin(phi)};
   HiSpID_Point out;if(HiSpID_sample(s,1,x,&out))return 3;
   J X=J::variable(x[0],1),Y=J::variable(x[1],2),Z=J::variable(x[2],3);
   J rp=hispid::sqrt((X-J(3))*(X-J(3))+Y*Y+Z*Z),rm=hispid::sqrt((X+J(3))*(X+J(3))+Y*Y+Z*Z);
   J cx=(rp+rm)/J(6),e=(rm-rp)/J(6),q=(cx-J(1))/(cx+J(1));
   J za=J(2)*q/(J(lambda)+J(1-lambda)*q)-J(1),arg=e*J(std::tanh((long double)kappa));
   J ze=hispid::unary(arg,std::atanh(arg.v),1/(1-arg.v*arg.v),2*arg.v/std::pow(1-arg.v*arg.v,2))/J(kappa);
   J P=J(100)+J(1e-6)*cheb(na-1,za)*cheb(nb-1,ze);
   J poly=Y*Y*Y*Y-J(6)*Y*Y*Z*Z+Z*Z*Z*Z;
   J beta=-hispid::power(J(1)-q,5)*poly*P*J(std::sqrt(2./16)/(8*81));
   hispid::Background bg{};if(hispid::background_geometry(c,x,bg,true,HISPID_SEED_TRUMPET_R0_M))return 4;
   hispid::CachedFields f{};f[2][0]=beta.v;for(int d=0;d<3;d++)f[2][d+1]=beta.d[d+1];
   double L[3][3];hispid::L_and_div(bg.metric,bg.inv,bg.C,f,L,nullptr);
   value_error=std::max(value_error,std::abs(out.correction[2]-double(beta.v))/(1+std::abs(double(beta.v))));
   for(int i=0;i<3;i++)for(int j=0;j<3;j++){
    const double psi=bg.psi.v,expected=(double(bg.M[i][j].v)+L[i][j])/(psi*psi)+std::pow(psi,4)*double(bg.metric[i][j].v*bg.K.v)/3;
    K_error=std::max(K_error,std::abs(out.Kij[3*i+j]-expected)/(1+std::abs(expected)));
   }
  }
  HiSpID_destroy(s);
 }
 std::cout<<"highest-degree perturbed constant: correction_error="<<value_error<<" K_error="<<K_error<<'\n';
 return value_error<1e-10&&K_error<1e-9?0:5;
}
