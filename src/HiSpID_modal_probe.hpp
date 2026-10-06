#ifndef HISPID_MODAL_PROBE_HPP
#define HISPID_MODAL_PROBE_HPP
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <cmath>
namespace hispid {
// Opt-in diagnostic only: preserve the assembled, unfactored sparse matrix
// and test its inverse with an all-frequency deterministic witness.
struct ModalFactorProbe {
 std::vector<int>offset,column;
 std::vector<double>value,exact,rhs;
 bool enabled=false;
 void capture(int na,int nb,const std::vector<double>&diagonal,
              const std::vector<double>&lower,const std::vector<double>&upper,
              const std::vector<double>&endpoint){
  const char*setting=std::getenv("HISPID_PROBE_MODAL_FACTORS");
  enabled=setting&&std::strcmp(setting,"1")==0;if(!enabled)return;
  const int n=na*nb;exact.resize(n);rhs.resize(n);offset.push_back(0);
  for(int row=0;row<n;row++)exact[row]=std::sin(1.71*row)+.2*std::cos(.31*row);
  auto add=[&](int c,double v){if(v!=0){column.push_back(c);value.push_back(v);}};
  for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
   const int row=j*na+i;
   if(!endpoint.empty()&&(j==0||j==nb-1)){
    for(int k=0;k<nb;k++)add(k*na+i,endpoint[(j==0?0:nb)+k]);
   }else{
    for(int k=0;k<na;k++)add(j*na+k,diagonal[size_t(row)*na+k]);
    if(j)add(row-na,lower[row]);if(j+1<nb)add(row+na,upper[row]);
   }
   long double sum=0;for(size_t k=offset.back();k<value.size();k++)sum+=(long double)value[k]*exact[column[k]];
   rhs[row]=double(sum);offset.push_back(int(value.size()));
  }
 }
 void report(int group,const std::vector<double>&solution)const{
  if(!enabled)return;
  double forward=0,backward=0;long double residual2=0,rhs2=0;
  for(size_t row=0;row<exact.size();row++){
   forward=std::max(forward,std::abs(solution[row]-exact[row]));
   long double sum=0,scale=std::abs(rhs[row]);
   for(int k=offset[row];k<offset[row+1];k++){
    long double term=(long double)value[k]*solution[column[k]];sum+=term;scale+=std::abs(term);
   }
   long double error=sum-rhs[row];residual2+=error*error;rhs2+=(long double)rhs[row]*rhs[row];
   backward=std::max(backward,double(std::abs(error)/(scale+1e-300L)));
  }
  std::fprintf(stderr,"HiSpID modal probe group=%d n=%zu forward_linf=%.17g backward_max=%.17g relative_l2=%.17g\n",
    group,exact.size(),forward,backward,double(std::sqrt(residual2/rhs2)));
  std::fflush(stderr);
 }
};
}
#endif
