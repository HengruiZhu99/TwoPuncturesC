#ifndef HISPID_RADIAL_PRECONDITIONER_HPP
#define HISPID_RADIAL_PRECONDITIONER_HPP
#include "HiSpID_spectral.hpp"
namespace hispid {
// Exact polynomial radial differentiation, used only in the approximate
// inverse. Dense radial blocks already exist; polar coupling stays local.
struct RadialPreconditioner {
 int n=0;
 std::vector<double>first,second,t,tz,tzz;
 explicit RadialPreconditioner(int count,double lambda):n(count),first(n*n),second(n*n),t(n),tz(n),tzz(n){
  for(int i=0;i<n;i++){
   const double z=-std::cos(std::acos(-1.)*(i+.5)/n),sigma=.5*(1+z);
   const double den=1-(1-lambda)*sigma;
   t[i]=lambda*sigma/den;tz[i]=.5*lambda/(den*den);
   tzz[i]=.5*lambda*(1-lambda)/(den*den*den);
   long double sum1=0,sum2=0;
   for(int j=0;j<n;j++)if(i!=j){
    const auto d=chebyshev_offdiagonal<long double>(n,i,j,1);
    const auto dd=chebyshev_offdiagonal<long double>(n,i,j,2);
    first[i*n+j]=double(d);second[i*n+j]=double(dd);sum1+=d;sum2+=dd;
   }
   first[i*n+i]=-double(sum1);second[i*n+i]=-double(sum2);
  }
 }
 double entry(int i,int j,int exponent)const{
  const double a=t[i],d=tz[i],dd=tzz[i];
  const double ctt=a*(1-a)*(1-a),ct=(1-a)*((exponent+1)*(1-a)-2*a);
  return ctt/(d*d)*second[i*n+j]+(ct/d-ctt*dd/(d*d*d))*first[i*n+j];
 }
};
}
#endif
