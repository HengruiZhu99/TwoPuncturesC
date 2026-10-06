#ifndef PUNCTURE_MODAL_PROJECTION_HPP
#define PUNCTURE_MODAL_PROJECTION_HPP
#include <cmath>
#ifdef PUNCTURES_KOKKOS
#include <Kokkos_Core.hpp>
#define PUNCTURE_PROJECTION_INLINE KOKKOS_INLINE_FUNCTION
#else
#define PUNCTURE_PROJECTION_INLINE inline
#endif
namespace puncture {
PUNCTURE_PROJECTION_INLINE void projection_two_sum(double a,double b,double&s,double&e){
 s=a+b;double q=s-a;e=(a-(s-q))+(b-q);
}
// Compensate subtraction, multiplication and summation. Tiny Fourier modes
// are divided by near-axis regularity factors; ordinary dot-product roundoff
// can consequently become a large preconditioned correction.
PUNCTURE_PROJECTION_INLINE double compensated_projection(
 const double*x,int stride,const double*w,int n,double mean){
 double sum=0,correction=0;
 for(int k=0;k<n;k++){
  double high,low;projection_two_sum(x[k*stride],-mean,high,low);
  volatile double rounded_product=w[k]*high;
  double product=rounded_product,error=::fma(w[k],high,-product)+w[k]*low;
  double next,addition_error;projection_two_sum(sum,product,next,addition_error);
  correction+=addition_error+error;sum=next;
 }
 return sum+correction;
}
}
#undef PUNCTURE_PROJECTION_INLINE
#endif
