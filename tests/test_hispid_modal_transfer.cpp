#include "../src/HiSpID_modal_transfer.hpp"
#include <chrono>
#include <cmath>
#include <iostream>
int main(){
 bool passed=true;
 for(size_t n:{size_t(64),size_t(512)}){
  std::vector<double>a(n*n),lu,upper(n),columns(n*n),batch(n*n),scratch;
  // Permute a diagonally dominant dense matrix to force nontrivial LU pivots.
  for(size_t i=0;i<n;i++)for(size_t j=0;j<n;j++){
   size_t row=(i+7)%n;
   a[i*n+j]=(row==j?2.:0.)+.1*std::sin(double(3*row+5*j+1))/n;
  }
  for(size_t i=0;i<n;i++)upper[i]=.3+.2*std::cos(double(i));
  lu=a;auto A=gsl_matrix_view_array(lu.data(),n,n);
  auto*p=gsl_permutation_alloc(n);int sign;gsl_linalg_LU_decomp(&A.matrix,p,&sign);
  auto now=[](){return std::chrono::steady_clock::now();};
  auto start=now();hispid::modal_transfer(&A.matrix,p,upper.data(),columns.data(),false,scratch);
  double old_seconds=std::chrono::duration<double>(now()-start).count();
  start=now();hispid::modal_transfer(&A.matrix,p,upper.data(),batch.data(),true,scratch);
  double new_seconds=std::chrono::duration<double>(now()-start).count();
  double difference=0,residual=0;bool pivoted=false;
  for(double x:columns)passed &= std::isfinite(x);
  for(double x:batch)passed &= std::isfinite(x);
  for(size_t i=0;i<n;i++){
   pivoted|=p->data[i]!=i;
   for(size_t j=0;j<n;j++){
    difference=std::max(difference,std::abs(batch[i*n+j]-columns[i*n+j])/(1+std::abs(columns[i*n+j])));
    long double sum=0,scale=0;
    for(size_t k=0;k<n;k++){long double term=(long double)a[i*n+k]*batch[k*n+j];sum+=term;scale+=std::abs(term);}
    double rhs=i==j?upper[i]:0;
    residual=std::max(residual,double(std::abs(sum-rhs)/(1+scale+std::abs(rhs))));
   }
  }
  passed &= pivoted&&std::isfinite(difference)&&std::isfinite(residual)&&difference<1e-12&&residual<1e-12;
  std::cout<<n<<" "<<difference<<" "<<residual<<" "<<old_seconds<<" "<<new_seconds<<" "<<pivoted<<"\n";
  gsl_permutation_free(p);
 }
 return passed?0:1;
}
