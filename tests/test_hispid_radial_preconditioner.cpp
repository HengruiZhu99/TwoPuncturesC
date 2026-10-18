#include "../src/HiSpID_radial_preconditioner.hpp"
#include <iostream>
int main(){
 double worst=0;
 for(int n:{12,32})for(double lambda:{.03,.2,1.})for(int r:{0,1,4}){
  hispid::RadialPreconditioner op(n,lambda);
  for(int degree:{0,3,n-1})for(int i=0;i<n;i++){
   const long double pi=std::acos(-1.L),z=-std::cos(pi*(i+.5L)/n);
   const long double theta=std::acos(z),u=std::cos(degree*theta);
   const long double du=degree*std::sin(degree*theta)/std::sin(theta);
   const long double ddu=(z*du-degree*degree*u)/(1-z*z);
   const long double s=(1+z)/2,d=1-(1-lambda)*s;
   long double t=lambda*s/d,tz=lambda/(2*d*d),tzz=lambda*(1-lambda)/(2*d*d*d);
   if constexpr(hispid::AxisDerivatives::dual_endpoint_map){
    const long double b=(2-lambda)/5,c=1-lambda-b,q=1-s;
    t=lambda*s+b*std::pow(s,5)+c*std::pow(s,15)*(136*q*q+17*s*q+s*s);
    tz=(lambda+5*b*std::pow(s,4)+2040*c*std::pow(s,14)*q*q)/2;
    tzz=(20*b*std::pow(s,3)+2040*c*(14*std::pow(s,13)*q*q-2*std::pow(s,14)*q))/4;
   }else if constexpr(hispid::AxisDerivatives::radial_polynomial){
    t=lambda*s+(1-lambda)*s*s;tz=(lambda+2*(1-lambda)*s)/2;tzz=(1-lambda)/2;
   }else if constexpr(hispid::AxisDerivatives::infinity_power==2){
    const long double q=lambda*s+(1-s)*(1-s),dq=lambda-2+2*s;
    t=lambda*s/q;tz=lambda*(1-s*s)/(2*q*q);
    tzz=-lambda*s/(2*q*q)-lambda*(1-s*s)*dq/(2*q*q*q);
   }
   const long double ctt=t*(1-t)*(1-t),ct=(1-t)*((r+1)*(1-t)-2*t);
   const long double exact=ctt*(ddu/(tz*tz)-du*tzz/(tz*tz*tz))+ct*du/tz;
   long double got=0;
   for(int j=0;j<n;j++)if(j!=i){
    const long double zj=-std::cos(pi*(j+.5L)/n);
    got+=op.entry(i,j,r)*(std::cos(degree*std::acos(zj))-u);
   }
   worst=std::max(worst,double(std::abs(got-exact)/(1+std::abs(exact))));
  }
 }
 std::cout<<"mapped polynomial radial operator worst scaled error="<<worst<<'\n';
 return worst<1e-9?0:1;
}
