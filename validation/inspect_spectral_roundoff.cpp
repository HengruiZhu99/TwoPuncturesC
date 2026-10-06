// Focused diagnostic: existing D*D versus an analytic barycentric D2.
// No production code or acceptance threshold is changed by this program.
#include "../src/HiSpID_spectral.hpp"
#include <iomanip>
#include <iostream>
int main(){
 std::cout<<std::setprecision(17);
 for(int n:{192,256,384,512}){
  int shape[3]={n,4,4};hispid::SpectralDerivatives old;old.initialize(shape);
  std::vector<long double>x(n),w(n),D(n*n),DD(n*n);
  const long double pi=acosl(-1.L);
  for(int i=0;i<n;i++){long double t=pi*(i+.5L)/n;x[i]=-cosl(t);w[i]=(i%2?-1.L:1.L)*sinl(t);}
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j){D[i*n+j]=w[j]/(w[i]*(x[i]-x[j]));D[i*n+i]-=D[i*n+j];}
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j){DD[i*n+j]=2*D[i*n+j]*(D[i*n+i]-1/(x[i]-x[j]));DD[i*n+i]-=DD[i*n+j];}
  for(int degree:{2,17}){
   std::vector<double>f(n);std::vector<long double>exact(n);
   for(int i=0;i<n;i++){
    // Match actual native double nodes, not the higher-precision construction nodes.
    long double z=-std::cos(std::acos(-1.)*(i+.5)/n),v0=1,v1=z,d0=0,d1=1,h0=0,h1=0;
    for(int k=2;k<=degree;k++){long double v=2*z*v1-v0,d=2*v1+2*z*d1-d0,h=4*d1+2*z*h1-h0;v0=v1;v1=v;d0=d1;d1=d;h0=h1;h1=h;}
    f[i]=double(v1);exact[i]=h1;
   }
   double eold=0,enew=0,innerold=0,innernew=0;
   for(int i=0;i<n;i++){
    double a=0,b=0;
    for(int j=0;j<n;j++)if(i!=j){a+=old.D2[0][i*n+j]*(f[j]-f[i]);b+=double(DD[i*n+j])*(f[j]-f[i]);}
    double ea=double(fabsl(a-exact[i])/(1+fabsl(exact[i]))),eb=double(fabsl(b-exact[i])/(1+fabsl(exact[i])));
    eold=std::max(eold,ea);enew=std::max(enew,eb);
    if(i>=n/8&&i<7*n/8){innerold=std::max(innerold,ea);innernew=std::max(innernew,eb);}
   }
   std::cout<<n<<' '<<degree<<' '<<eold<<' '<<enew<<' '<<innerold<<' '<<innernew<<'\n';
  }
 }
}
