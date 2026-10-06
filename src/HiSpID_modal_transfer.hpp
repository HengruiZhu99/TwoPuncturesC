#ifndef HISPID_MODAL_TRANSFER_HPP
#define HISPID_MODAL_TRANSFER_HPP
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#ifndef HISPID_BATCHED_MODAL_TRANSFER
#define HISPID_BATCHED_MODAL_TRANSFER 0
#endif
namespace hispid {
// Form A^{-1} diag(upper) from PA=LU. The optional matrix-RHS solve avoids
// one strided triangular solve per column; both paths form the same transfer.
inline void modal_transfer(const gsl_matrix *lu,const gsl_permutation *p,
                           const double *upper,double *out,bool batched,
                           std::vector<double>&scratch){
 const size_t n=lu->size1;
 if(batched){
  std::fill(out,out+n*n,0.);
  for(size_t i=0;i<n;i++)out[i*n+p->data[i]]=upper[p->data[i]];
  auto rhs=gsl_matrix_view_array(out,n,n);
  if(gsl_blas_dtrsm(CblasLeft,CblasLower,CblasNoTrans,CblasUnit,1.,lu,&rhs.matrix)||
     gsl_blas_dtrsm(CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,1.,lu,&rhs.matrix))
   throw std::runtime_error("Batched modal transfer solve failed");
 }else{
  scratch.resize(n);
  for(size_t q=0;q<n;q++){
   std::fill(scratch.begin(),scratch.end(),0.);scratch[q]=upper[q];
   auto b=gsl_vector_view_array(scratch.data(),n);
   auto x=gsl_vector_view_array_with_stride(out+q,n,n);
   if(gsl_linalg_LU_solve(lu,p,&b.vector,&x.vector))
    throw std::runtime_error("Modal transfer solve failed");
  }
 }
}
}
#endif
