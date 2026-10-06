#include "../src/HiSpID_modal_block.hpp"
#include <iostream>
#include <iomanip>
// Independent endpoint evaluator via Chebyshev coefficients, not the
// barycentric weights used to construct the border correction.
static double endpoint(int n,int j,int side){
 double value=1./n,theta=std::acos(-1.)*(j+.5)/n;
 for(int k=1;k<n;k++)value+=2./n*std::cos(k*theta)*(side&&k%2?-1.:1.);
 return value;
}
static std::vector<double> barycentric(int n){
 std::vector<double>w(2*n);const double pi=std::acos(-1.);
 for(int s=0;s<2;s++){
  double sum=0;
  for(int j=0;j<n;j++){double t=pi*(j+.5)/n;w[s*n+j]=(j%2?-1.:1.)*std::sin(t)/((s?1.:-1.)+std::cos(t));sum+=w[s*n+j];}
  for(int j=0;j<n;j++)w[s*n+j]/=sum;
 }return w;
}
int main(){
 gsl_set_error_handler_off();bool passed=true;
 for(auto dims:{std::pair<int,int>{5,8},{12,24}}){
  int a=dims.first,b=dims.second,N=a*b;
  hispid::ModalBlock B;B.na=a;B.nb=b;
  B.lu.assign(a*a*b,0.);B.transfer.resize(a*a*b);B.lower.assign(N,0.);B.upper.assign(N,0.);B.permutation.resize(N);
  std::vector<double>A(size_t(N)*N,0.);
  for(int j=0;j<b;j++)for(int i=0;i<a;i++){
   int row=j*a+i;
   if(j==0||j==b-1){
    int next=j==0?1:b-2;double t=std::acos(-1.)/b,z0=-std::cos(t*(j+.5)),z1=-std::cos(t*(next+.5)),end=j==0?-1.:1.;
    B.lu[(j*a+i)*a+i]=(z1-end)/(z1-z0);
    (j==0?B.upper:B.lower)[row]=(end-z0)/(z1-z0);
   }else if(i==0){
    for(int k=0;k<a;k++)B.lu[(j*a+i)*a+k]=endpoint(a,k,0);
   }else{
    B.lu[(j*a+i)*a+i]=-9.-.1*i;
    B.lu[(j*a+i)*a+i-1]=3.;
    if(i+1<a)B.lu[(j*a+i)*a+i+1]=1.2;
    B.lower[row]=.8;B.upper[row]=1.1;
   }
   for(int k=0;k<a;k++)A[size_t(row)*N+j*a+k]=B.lu[(j*a+i)*a+k];
   if(j)A[size_t(row)*N+row-a]=B.lower[row];
   if(j+1<b)A[size_t(row)*N+row+a]=B.upper[row];
  }
  B.factor();bool pivoted=false;
  for(int j=0;j<b;j++)for(int i=0;i<a;i++)pivoted|=B.permutation[j*a+i]!=size_t(i);
  // Dense target differs in its polar rows, retaining radial corner priority.
  for(int side=0;side<2;side++)for(int i=0;i<a;i++){
   int row=(side?b-1:0)*a+i;
   std::fill_n(A.data()+size_t(row)*N,N,0.);
   for(int j=0;j<b;j++)A[size_t(row)*N+j*a+i]=endpoint(b,j,side);
  }
  B.prepare_polar(barycentric(b));
  std::vector<double>exact(N),rhs(N),solution(4*N*2,0.);
  // Nonzero inhomogeneous boundary data and all nodal frequencies: a zero
  // boundary polynomial alone would miss errors in U and the Schur solve.
  for(int i=0;i<N;i++)exact[i]=std::sin(1.71*i)+.2*std::cos(.31*i);
  for(int i=0;i<N;i++)for(int j=0;j<N;j++)rhs[i]+=A[size_t(i)*N+j]*exact[j];
  for(int i=0;i<N;i++)solution[4*(i+N)+2]=rhs[i];
  B.solve(solution,1,2);
  double error=0,residual=0,leak=0;
  for(int i=0;i<N;i++){
   double x=solution[4*(i+N)+2];passed&=std::isfinite(x);error=std::max(error,std::abs(x-exact[i]));
   long double sum=0,scale=1+std::abs(rhs[i]);
   for(int j=0;j<N;j++){long double term=(long double)A[size_t(i)*N+j]*solution[4*(j+N)+2];sum+=term;scale+=std::abs(term);}
   residual=std::max(residual,double(std::abs(sum-rhs[i])/scale));
  }
  for(int i=0;i<8*N;i++)if(i<4*N||i%4!=2)leak=std::max(leak,std::abs(solution[i]));
  passed&=pivoted&&error<1e-11&&residual<1e-12&&leak==0;
  std::cout<<a<<"x"<<b<<" pivoted="<<pivoted<<" solution_error="<<std::setprecision(15)<<error<<" scaled_residual="<<residual<<" interleaved_leak="<<leak<<"\n";
 }
 return passed?0:1;
}
