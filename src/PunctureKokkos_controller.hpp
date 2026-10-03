/* Device vector policy for the shared C-compatible Krylov controller. Only
 * vector banks are device allocations; every scalar/pointer bank stays host. */
#ifndef PUNCTURE_KOKKOS_CONTROLLER_HPP
#define PUNCTURE_KOKKOS_CONTROLLER_HPP
#include "PunctureKokkos.hpp"
#include <cstdint>
#include <climits>
#include <unordered_map>
namespace puncture { namespace controller {
using std::isfinite;using std::sqrt;using std::hypot;using std::fabs;
static size_t live_workspace=0;
static std::unordered_map<double*,size_t>banks;
static void *bank(size_t rows,size_t cols,size_t width){if(cols&&rows>SIZE_MAX/cols)return nullptr;size_t n=rows*cols;if(width&&n>SIZE_MAX/width)return nullptr;return calloc(n,width);}
static double*vector_bank(size_t n){
 if(n>SIZE_MAX/sizeof(double))return nullptr;
 double*v=nullptr;
 try{v=static_cast<double*>(Kokkos::kokkos_malloc<typename Exec::memory_space>("Krylov workspace",n*sizeof(double)));if(!v)return nullptr;banks.emplace(v,n*sizeof(double));live_workspace+=n*sizeof(double);
  statistics().peak_workspace_bytes=std::max(statistics().peak_workspace_bytes,(unsigned long long)live_workspace);return v;
 }catch(...){if(v)Kokkos::kokkos_free<typename Exec::memory_space>(v);return nullptr;}
}
static void vector_free(double*v){if(v){Exec().fence();auto it=banks.find(v);if(it!=banks.end()){live_workspace-=it->second;banks.erase(it);}Kokkos::kokkos_free<typename Exec::memory_space>(v);}}
static void vector_zero(double*v,int n){Kokkos::parallel_for("zero",Range(0,n),KOKKOS_LAMBDA(int i){v[i]=0;});}
static void vector_difference(double*r,const double*b,const double*Ax,int n){Kokkos::parallel_for("difference",Range(0,n),KOKKOS_LAMBDA(int i){r[i]=b[i]-Ax[i];});}
static void vector_divide(double*out,const double*in,double beta,int n){Kokkos::parallel_for("divide",Range(0,n),KOKKOS_LAMBDA(int i){out[i]=in[i]/beta;});}
static void vector_subtract(double*w,const double*v,double a,int n){Kokkos::parallel_for("MGS update",Range(0,n),KOKKOS_LAMBDA(int q){w[q]-=a*v[q];});}
static void vector_add(double*x,const double*v,double a,int n){Kokkos::parallel_for("GMRES update",Range(0,n),KOKKOS_LAMBDA(int q){x[q]+=v[q]*a;});}
static void vector_add_left(double*x,const double*v,double a,int n){Kokkos::parallel_for("BiCG early update",Range(0,n),KOKKOS_LAMBDA(int j){x[j]+=a*v[j];});}
static void vector_shadow(double*rt,double*r,const double*b,int n){Kokkos::parallel_for("shadow residual",Range(0,n),KOKKOS_LAMBDA(int j){rt[j]=r[j]=b[j]-r[j];});}
static void vector_copy(double*p,const double*r,int n){Kokkos::parallel_for("copy",Range(0,n),KOKKOS_LAMBDA(int j){p[j]=r[j];});}
static void vector_search(double*p,const double*r,const double*v,double beta,double omega,int n){Kokkos::parallel_for("search vector",Range(0,n),KOKKOS_LAMBDA(int j){p[j]=r[j]+beta*(p[j]-omega*v[j]);});}
static void vector_step(double*s,const double*r,const double*v,double alpha,int n){Kokkos::parallel_for("BiCG half residual",Range(0,n),KOKKOS_LAMBDA(int j){s[j]=r[j]-alpha*v[j];});}
static void vector_bicg_update(double*x,double*r,const double*s,const double*t,const double*ph,const double*sh,double alpha,double omega,int n){Kokkos::parallel_for("BiCG update",Range(0,n),KOKKOS_LAMBDA(int j){x[j]+=alpha*ph[j]+omega*sh[j];r[j]=s[j]-omega*t[j];});}
static double dot(const double*a,const double*b,int n){double q=0;Kokkos::parallel_reduce("Krylov dot",Range(0,n),KOKKOS_LAMBDA(int i,double&v){v+=a[i]*b[i];},q);return q;}
static double length(const double*a,int n){return sqrt(dot(a,a,n));}
static double inner(const PK_Options*,const double*a,const double*b,int n){return dot(a,b,n);}
static double norm(const PK_Options*,const double*a,int n){return length(a,n);}
static int finite_vector(const double*a,int n){int bad=0;Kokkos::parallel_reduce("Krylov finite",Range(0,n),KOKKOS_LAMBDA(int i,int&v){v+=!Kokkos::isfinite(a[i]);},bad);return bad==0;}
struct Context {const Action&A;const Action&M;int n;};
static int apply_A(void*c,const double*b,double*x){auto&s=*static_cast<Context*>(c);try{s.A(View(const_cast<double*>(b),s.n),View(x,s.n));return 0;}catch(...){return -1;}}
static int apply_M(void*c,const double*b,double*x){auto&s=*static_cast<Context*>(c);try{s.M(View(const_cast<double*>(b),s.n),View(x,s.n));return 0;}catch(...){return -1;}}
static int action(PK_Apply f,void*c,const double*b,double*x,int n,int*count){
 (*count)++;Timed t(f==apply_A?statistics().operator_seconds:statistics().precondition_seconds);if(f(c,b,x))return PK_CALLBACK;return finite_vector(x,n)?PK_SUCCESS:PK_NONFINITE;
}
#include "PunctureKrylov_controller.inc"
} }
#endif
