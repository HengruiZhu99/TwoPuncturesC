#ifndef HISPID_MODAL_BLOCK_HPP
#define HISPID_MODAL_BLOCK_HPP
#include "HiSpID_modal_transfer.hpp"
#include "HiSpID_modal_probe.hpp"
#include <cmath>
namespace hispid {
struct ModalBlock {
 int na=0,nb=0;
 std::vector<double>polar_endpoint;
 std::vector<double>lu,transfer,lower,upper;
 std::vector<size_t>permutation;
 // Optional exact polar endpoint correction, built after the base factors.
 std::vector<double>polar_delta,polar_response,polar_inverse;
 void factor(int diagnostic_group=-1){
  ModalFactorProbe probe;probe.capture(na,nb,lu,lower,upper,polar_endpoint);
  polar_delta.clear();polar_response.clear();polar_inverse.clear();
  std::vector<double>rhs(na);int sign=0;
  for(int j=0;j<nb;j++){
   double*block=lu.data()+j*na*na;
   if(j)for(int i=0;i<na;i++)for(int q=0;q<na;q++)
    block[i*na+q]-=lower[j*na+i]*transfer[((j-1)*na+i)*na+q];
   auto A=gsl_matrix_view_array(block,na,na);
   gsl_permutation p{size_t(na),permutation.data()+j*na};gsl_permutation_init(&p);
   if(gsl_linalg_LU_decomp(&A.matrix,&p,&sign))throw std::runtime_error("Modal block factorization failed");
   for(int i=0;i<na;i++)if(!std::isfinite(block[i*na+i])||std::abs(block[i*na+i])<1e-30)
    throw std::runtime_error("Modal block singular pivot");
   if(j+1<nb)hispid::modal_transfer(&A.matrix,&p,upper.data()+j*na,
       transfer.data()+j*na*na,bool(HISPID_BATCHED_MODAL_TRANSFER),rhs);
  }
  if(!polar_endpoint.empty())prepare_polar(polar_endpoint);
  if(probe.enabled){auto x=probe.rhs;solve_many(x.data(),1);correct_polar(x.data());probe.report(diagnostic_group,x);}
 }
 // Solve many contiguous row-major RHS columns with the existing factors.
 // Used to construct A0^{-1} U without one triangular solve per column.
 void solve_many(double*x,int columns)const{
  std::vector<double>permuted(size_t(na)*columns);
  for(int j=0;j<nb;j++){
   double*row=x+size_t(j)*na*columns;
   if(j)for(int i=0;i<na;i++)for(int k=0;k<columns;k++)
    row[i*columns+k]-=lower[j*na+i]*row[(i-na)*columns+k];
   for(int i=0;i<na;i++)std::copy_n(row+permutation[j*na+i]*columns,columns,
                                  permuted.data()+i*columns);
   std::copy(permuted.begin(),permuted.end(),row);
   auto A=gsl_matrix_const_view_array(lu.data()+size_t(j)*na*na,na,na);
   auto R=gsl_matrix_view_array(row,na,columns);
   if(gsl_blas_dtrsm(CblasLeft,CblasLower,CblasNoTrans,CblasUnit,1.,&A.matrix,&R.matrix)||
      gsl_blas_dtrsm(CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,1.,&A.matrix,&R.matrix))
    throw std::runtime_error("Modal matrix forward solve failed");
  }
  for(int j=nb-2;j>=0;j--){
   auto T=gsl_matrix_const_view_array(transfer.data()+size_t(j)*na*na,na,na);
   auto next=gsl_matrix_const_view_array(x+size_t(j+1)*na*columns,na,columns);
   auto row=gsl_matrix_view_array(x+size_t(j)*na*columns,na,columns);
   if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1.,&T.matrix,&next.matrix,1.,&row.matrix))
    throw std::runtime_error("Modal matrix backward solve failed");
  }
 }
 // Replace the two local polar extrapolation rows by exact endpoint rows.
 // A=A0+U V^T; response=A0^{-1}U; inverse=(I+V^T response)^{-1}.
 // endpoint contains the two nb-element endpoint interpolation rows.
 void prepare_polar(const std::vector<double>&endpoint){
  if(na<1||nb<2||endpoint.size()<size_t(2*nb))
   throw std::runtime_error("Invalid polar border dimensions");
  const int rank=2*na;const size_t N=size_t(na)*nb;
  polar_delta.assign(endpoint.begin(),endpoint.begin()+2*nb);
  for(int side=0;side<2;side++){
   int j=side?nb-1:0,next=side?nb-2:1;
   const double pi=std::acos(-1.),z0=-std::cos(pi*(j+.5)/nb),
    z1=-std::cos(pi*(next+.5)/nb),end=side?1.:-1.;
   polar_delta[side*nb+j]-=(z1-end)/(z1-z0);
   polar_delta[side*nb+next]-=(end-z0)/(z1-z0);
  }
  polar_response.assign(N*rank,0.);
  for(int side=0;side<2;side++)for(int i=0;i<na;i++)
   polar_response[(size_t(side?nb-1:0)*na+i)*rank+side*na+i]=1.;
  solve_many(polar_response.data(),rank);
  std::vector<double>schur(size_t(rank)*rank,0.);
  for(int side=0;side<2;side++)for(int i=0;i<na;i++)for(int k=0;k<rank;k++){
   double value=(side*na+i==k)?1.:0.;
   for(int j=0;j<nb;j++)value+=polar_delta[side*nb+j]*polar_response[(size_t(j)*na+i)*rank+k];
   schur[(side*na+i)*rank+k]=value;
  }
  auto S=gsl_matrix_view_array(schur.data(),rank,rank);
  std::vector<size_t>pivots(rank);gsl_permutation p{size_t(rank),pivots.data()};
  gsl_permutation_init(&p);int sign=0;
  if(gsl_linalg_LU_decomp(&S.matrix,&p,&sign))throw std::runtime_error("Polar Schur factorization failed");
  for(int i=0;i<rank;i++)if(!std::isfinite(schur[i*rank+i])||std::abs(schur[i*rank+i])<1e-30)
   throw std::runtime_error("Polar Schur singular pivot");
  polar_inverse.resize(size_t(rank)*rank);
  auto inv=gsl_matrix_view_array(polar_inverse.data(),rank,rank);
  if(gsl_linalg_LU_invert(&S.matrix,&p,&inv.matrix))throw std::runtime_error("Polar Schur inverse failed");
 }
 void correct_polar(double*x)const{
  if(polar_inverse.empty())return;
  const int rank=2*na;std::vector<double>g(rank,0.),z(rank,0.);
  for(int side=0;side<2;side++)for(int i=0;i<na;i++)for(int j=0;j<nb;j++)
   g[side*na+i]+=polar_delta[side*nb+j]*x[j*na+i];
  for(int i=0;i<rank;i++)for(int k=0;k<rank;k++)z[i]+=polar_inverse[i*rank+k]*g[k];
  for(size_t row=0;row<size_t(na)*nb;row++){
   double correction=0;for(int k=0;k<rank;k++)correction+=polar_response[row*rank+k]*z[k];
   x[row]-=correction;
  }
 }
 void solve(double*x,int mode,int component)const{
  std::vector<double>f(na*nb);
  for(int j=0;j<nb;j++){
   for(int i=0;i<na;i++)f[j*na+i]=x[4*(i+na*(j+nb*mode))+component]-(j?lower[j*na+i]*f[(j-1)*na+i]:0);
   auto A=gsl_matrix_const_view_array(lu.data()+j*na*na,na,na);
   gsl_permutation p{size_t(na),const_cast<size_t*>(permutation.data()+j*na)};
   auto y=gsl_vector_view_array(f.data()+j*na,na);
   if(gsl_linalg_LU_svx(&A.matrix,&p,&y.vector))throw std::runtime_error("Modal forward solve failed");
  }
  for(int j=nb-2;j>=0;j--)for(int i=0;i<na;i++){
   double sum=0;for(int q=0;q<na;q++)sum+=transfer[(j*na+i)*na+q]*f[(j+1)*na+q];f[j*na+i]-=sum;
  }
  correct_polar(f.data());
  for(int j=0;j<nb;j++)for(int i=0;i<na;i++)x[4*(i+na*(j+nb*mode))+component]=f[j*na+i];
 }
 void solve(std::vector<double>&x,int mode,int component)const{solve(x.data(),mode,component);}
};
}
#endif
