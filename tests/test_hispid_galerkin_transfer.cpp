#include "HiSpID_galerkin.hpp"
#include <iostream>
int main(){
 const double pi=std::acos(-1.);double roundtrip=0,interpolation=0,alias=0;
 for(int a:{8,16,31}){
  const int n[3]={a,2*a,4},q[3]={(3*a+1)/2,3*a,4};hispid::GalerkinTransfer t(n,q,1);
  std::vector<double>x(size_t(n[0])*n[1]*n[2]),fine(size_t(q[0])*q[1]*q[2]),back(x.size());
  for(int k=0;k<4;k++)for(int j=0;j<n[1];j++)for(int i=0;i<n[0];i++)x[i+n[0]*(j+n[1]*k)]=3.+(k+1)*std::cos((n[0]-1)*pi*(i+.5)/n[0])*std::cos((n[1]-1)*pi*(j+.5)/n[1]);
  t.prolong(x.data(),fine.data());t.restrict(fine.data(),back.data());
  for(size_t p=0;p<x.size();p++)roundtrip=std::max(roundtrip,std::abs(x[p]-back[p]));
  for(int k=0;k<4;k++)for(int j=0;j<q[1];j++)for(int i=0;i<q[0];i++){
   double exact=3.+(k+1)*std::cos((n[0]-1)*pi*(i+.5)/q[0])*std::cos((n[1]-1)*pi*(j+.5)/q[1]);
   interpolation=std::max(interpolation,std::abs(fine[i+q[0]*(j+q[1]*k)]-exact));
  }
  for(int k=0;k<4;k++)for(int j=0;j<q[1];j++)for(int i=0;i<q[0];i++)fine[i+q[0]*(j+q[1]*k)]=std::pow(std::cos((n[0]-1)*pi*(i+.5)/q[0]),2);
  t.restrict(fine.data(),back.data());for(double z:back)alias=std::max(alias,std::abs(z-.5));
 }
 std::cout<<"interpolation="<<interpolation<<" roundtrip="<<roundtrip<<" quadratic_alias="<<alias<<'\n';
 return interpolation<5e-12&&roundtrip<5e-12&&alias<5e-12?0:1;
}
