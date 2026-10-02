// Independent automatic-differentiation control of the factored modal basis.
#include <iostream>
#include "../src/HiSpID_axis.hpp"
#include "../src/HiSpID_jets.hpp"
int main(){
 double error=0,orthogonal_error=0;const double pi=std::acos(-1.0);
 const double lambda=hispid::AxisDerivatives::radial_stretch,kappa=hispid::AxisDerivatives::angular_stretch;
 for(int np:{8,16})for(int m=0;m<=np/2;m++)for(bool sine:{false,true}){
  if(sine&&(m==0||m==np/2))continue;
  const int shape[3]={12,13,np},count=shape[0]*shape[1]*np,mode=sine?np/2+m:m,r=hispid::AxisDerivatives::exponent(m);
  hispid::AxisDerivatives derivative;derivative.initialize(shape);
  std::array<std::vector<double>,10> storage;for(auto&v:storage)v.resize(count);
  derivs w{};double**targets[10]={&w.d0,&w.d1,&w.d2,&w.d3,&w.d11,&w.d12,&w.d13,&w.d22,&w.d23,&w.d33};
  for(int d=0;d<10;d++)*targets[d]=storage[d].data();
  const double normal=std::sqrt((m==0||m==np/2?1.:2.)/np);
  for(int j=0;j<shape[1];j++)for(int i=0;i<shape[0];i++){
   const double a=.5*(derivative.coordinate[0][i]+1),B=derivative.coordinate[1][j],t=a*a,eta=-2*B/(1+B*B),sigma=t/(lambda+(1-lambda)*t),zeta=std::atanh(eta*std::tanh(kappa))/kappa;
   storage[0][i+shape[0]*(j+shape[1]*mode)]=(1+.2*sigma+.07*zeta+.03*sigma*zeta+.02*sigma*sigma)/normal;
  }
  derivative.apply(1,&w);
  for(int k=0;k<np;k++)for(int j=0;j<shape[1];j++)for(int i=0;i<shape[0];i++){
   using hispid::Jet;Jet A=Jet::variable(derivative.coordinate[0][i],1),B=Jet::variable(derivative.coordinate[1][j],2),phi=Jet::variable(2*pi*k/np,3);
   Jet a=(A+Jet(1))/Jet(2),t=a*a,eta=Jet(-2)*B/(Jet(1)+B*B),sn=(Jet(1)-B*B)/(Jet(1)+B*B),factor(1);
   for(int exponent=0;exponent<r;exponent++)factor=factor*a*sn;
   Jet sigma=t/(Jet(lambda)+Jet(1-lambda)*t),angle=eta*Jet(std::tanh(kappa));Jet zeta=hispid::unary(angle,std::atanh(angle.v),1/(1-angle.v*angle.v),2*angle.v/std::pow(1-angle.v*angle.v,2))/Jet(kappa);
   Jet phase=Jet(m)*phi;const long double ss=std::sin(phase.v),cc=std::cos(phase.v);Jet wave=sine?hispid::unary(phase,ss,cc,-ss):hispid::unary(phase,cc,-ss,-cc);
   Jet V=(Jet(1)+a)*factor*(Jet(1)+Jet(.2)*sigma+Jet(.07)*zeta+Jet(.03)*sigma*zeta+Jet(.02)*sigma*sigma)*wave;
   const double expected[10]={double(V.v),double(V.d[1]),double(V.d[2]),double(V.d[3]),double(V.h[1][1]),double(V.h[1][2]),double(V.h[1][3]),double(V.h[2][2]),double(V.h[2][3]),double(V.h[3][3])};
   const int p=i+shape[0]*(j+shape[1]*k);
   for(int d=0;d<10;d++)error=std::max(error,std::abs(storage[d][p]-expected[d])/(1+std::abs(expected[d])));
  }
  for(int i=0;i<np;i++)for(int j=0;j<np;j++){
   double dot=0;for(int k=0;k<np;k++)dot+=derivative.forward[i*np+k]*derivative.forward[j*np+k];
   orthogonal_error=std::max(orthogonal_error,std::abs(dot-(i==j?1.:0.)));
  }
 }
 std::cout<<"factored modes versus independent AD (values/all nine derivatives) "<<error<<", Fourier orthogonality "<<orthogonal_error<<'\n';
 return error<2e-10&&orthogonal_error<2e-14?0:1;
}
