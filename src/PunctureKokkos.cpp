#include "PunctureKokkos.hpp"
#include "PunctureKokkos_controller.hpp"
#include "HiSpID_spectral.hpp"
#include "TP_Modal.h"
#include <gsl/gsl_complex_math.h>
#include <mutex>
#include <impl/Kokkos_Profiling.hpp>
namespace puncture {
namespace {
struct Runtime {
 bool owned=false;int threads=0;PunctureExecutionStats stats{};
 std::mutex memory_mutex;
 std::unordered_map<const void*,std::pair<uint64_t,bool>>allocations;
 Kokkos::Tools::Experimental::EventSet prior_callbacks{};
};
Runtime&runtime(){static Runtime r;return r;}
void allocation(Kokkos::Tools::SpaceHandle space,const char*label,const void*ptr,uint64_t bytes){
 auto&r=runtime();bool pinned=std::strstr(space.name,"Pinned")!=nullptr;
 bool device=!pinned&&(std::strstr(space.name,"Cuda")||std::strstr(space.name,"HIP")||std::strstr(space.name,"SYCL"));
 {std::lock_guard<std::mutex>lock(r.memory_mutex);r.allocations[ptr]={bytes,device};auto&current=device?r.stats.kokkos_device_current_bytes:r.stats.kokkos_host_current_bytes;auto&peak=device?r.stats.kokkos_device_peak_bytes:r.stats.kokkos_host_peak_bytes;current+=bytes;peak=std::max(peak,current);}
 if(r.prior_callbacks.allocate_data)r.prior_callbacks.allocate_data(space,label,ptr,bytes);
}
void deallocation(Kokkos::Tools::SpaceHandle space,const char*label,const void*ptr,uint64_t bytes){
 auto&r=runtime();{std::lock_guard<std::mutex>lock(r.memory_mutex);auto it=r.allocations.find(ptr);if(it!=r.allocations.end()){auto&current=it->second.second?r.stats.kokkos_device_current_bytes:r.stats.kokkos_host_current_bytes;current-=it->second.first;r.allocations.erase(it);}}
 if(r.prior_callbacks.deallocate_data)r.prior_callbacks.deallocate_data(space,label,ptr,bytes);
}
}
PunctureExecutionStats&statistics(){return runtime().stats;}
std::recursive_mutex&execution_mutex(){static std::recursive_mutex mutex;return mutex;}
void initialize(){if(Puncture_execution_initialize(0))throw std::runtime_error("Kokkos runtime unavailable or finalized");}
int solve(View b,View x,const PK_Options&o,const Action&A,const Action&M,PK_Result&result,PK_Monitor monitor,void*mc){
 ExecutionLock lock(execution_mutex());
 if(b.extent(0)>INT_MAX||b.extent(0)!=x.extent(0)||o.dot||o.norm){result={};return result.status=PK_INVALID;}
 Timed timer(statistics().linear_seconds);controller::Context context{A,M,int(b.extent(0))};
 // Legacy norm callbacks reference host memory and cannot enter this policy.
 PK_Options options=o;options.dot=nullptr;options.norm=nullptr;
 return controller::controller_solve(context.n,b.data(),x.data(),&options,A?controller::apply_A:nullptr,M?controller::apply_M:nullptr,&context,monitor,mc,&result);
}
void initialize_spectral(Operator&o,const int*n,double b,const std::vector<double>*coordinates){
 hispid::SpectralDerivatives raw;raw.initialize(n);
 for(int axis=0;axis<3;axis++){o.spectral.D[axis]=upload(raw.D[axis],"D");o.spectral.D2[axis]=upload(raw.D2[axis],"D2");}
 std::vector<double> coord[2];for(int axis=0;axis<2;axis++){
  if(coordinates)coord[axis]=coordinates[axis];else{coord[axis].resize(n[axis]);for(int i=0;i<n[axis];i++)coord[axis][i]=-std::cos(Pi*(i+.5)/n[axis]);}
  o.spectral.coordinate[axis]=upload(coord[axis],"coordinates");
 }
 std::vector<double>chain(size_t(n[0])*n[1]*6),trig(2*n[2]);
 for(int j=0;j<n[1];j++)for(int i=0;i<n[0];i++){
  double a=.5*(coord[0][i]+1),X=2*std::atanh(a),R=Pih+2*std::atan(coord[1][j]);
  auto C=gsl_complex_rect(X,R),c=gsl_complex_mul_real(gsl_complex_cosh(C),b),inv=gsl_complex_inverse(gsl_complex_mul_real(gsl_complex_sinh(C),b));
  auto cc=gsl_complex_negative(gsl_complex_mul(gsl_complex_mul(inv,inv),gsl_complex_mul(inv,c)));
  auto*q=chain.data()+6*(i+n[0]*j);q[0]=GSL_REAL(inv);q[1]=GSL_IMAG(inv);q[2]=GSL_REAL(cc);q[3]=GSL_IMAG(cc);q[4]=gsl_complex_abs2(inv);q[5]=1/GSL_IMAG(c);
 }
 for(int k=0;k<n[2];k++){trig[2*k]=std::cos(2*Pi*k/n[2]);trig[2*k+1]=std::sin(2*Pi*k/n[2]);}
 o.chain=upload(chain,"coordinate chain");o.trig=upload(trig,"azimuthal coordinates");
}
void Modal::prepare(){
#ifdef KOKKOS_ENABLE_CUDA
 if constexpr(std::is_same_v<Exec,Kokkos::Cuda>){
  const int a=na;auto L=lu;auto p=permutation;
  inverse=View(Kokkos::view_alloc(Kokkos::WithoutInitializing,"inverse LU"),lu.extent(0));auto inv=inverse;
  // Each column is independent; retain pivoting of the already factored CPU
  // Schur blocks. Column-major inverse makes CUDA row GEMV coalesced.
  Kokkos::parallel_for("invert radial blocks",Range(0,groups*nb*na),KOKKOS_LAMBDA(int task){
   int block=task/a,col=task%a,off=block*a*a,poff=block*a;
   for(int i=0;i<a;i++){double value=p(poff+i)==col?1:0;for(int q=0;q<i;q++)value-=L(off+i*a+q)*inv(off+col*a+q);inv(off+col*a+i)=value;}
   for(int i=a-1;i>=0;i--){double value=inv(off+col*a+i);for(int q=i+1;q<a;q++)value-=L(off+i*a+q)*inv(off+col*a+q);inv(off+col*a+i)=value/L(off+i*a+i);}
  });
  Exec().fence();if(!finite(inverse))throw std::runtime_error("nonfinite modal inverse");
  auto T=transfer,lower_=lower;const bool diagonal=diagonal_lower;
  // Transpose the existing banks in place. Neighboring CUDA row lanes now
  // load neighboring addresses, with no second bank or changed dot order.
  Kokkos::parallel_for("coalesce modal blocks",Range(0,groups*nb*na*na),KOKKOS_LAMBDA(int task){
   int block=task/(a*a),row=(task/a)%a,col=task%a,off=block*a*a;
   if(row<col){double value=T(off+row*a+col);T(off+row*a+col)=T(off+col*a+row);T(off+col*a+row)=value;
    if(!diagonal){value=lower_(off+row*a+col);lower_(off+row*a+col)=lower_(off+col*a+row);lower_(off+col*a+row)=value;}}
  });
  Exec().fence();
  lu=View();permutation=Indices();
 }
#endif
}
void Modal::apply(View in,View out,bool projected){
 const int a=na,b=nb,N=np,V=nv,stride=a*b,total=stride*N*V;
 auto f=forward,scale=row_scale,mod=column,y=workspace,L=lu,T=transfer,lower_=lower,inv=inverse;auto rows=block,perm=permutation;auto native=native_blocks;
 const bool diag=diagonal_lower,output_modal=modal_output;
 if(projected)copy(mod,in);
 else Kokkos::parallel_for("modal Fourier projection",Range(0,stride*V),KOKKOS_LAMBDA(int line){
  int v=line%V,row=line/V;double sum=0,error=0;
  for(int k=0;k<N;k++){double z=in(V*(row+stride*k)+v)-error,t=sum+z;error=(t-sum)-z;sum=t;}
  double mean=sum/N;
  for(int mode=0;mode<N;mode++){double value=mode==0?sum/Kokkos::sqrt(double(N)):0;
   if(mode)for(int k=0;k<N;k++)value+=f(mode*N+k)*(in(V*(row+stride*k)+v)-mean);
   int p=V*(row+stride*mode)+v;mod(p)=value/scale(p);
  }
 });
#ifdef KOKKOS_ENABLE_CUDA
 if constexpr(std::is_same_v<Exec,Kokkos::Cuda>){
  using Team=Kokkos::TeamPolicy<Exec>;using Member=Team::member_type;
  Kokkos::parallel_for("modal inverse CUDA",Team(N*V,128),KOKKOS_LAMBDA(const Member&team){
   int task=team.league_rank(),mode=task/V,v=task%V,g=rows(task);
   for(int j=0;j<b;j++){
    int off=(g*b+j)*a*a;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,a),[&](int i){int p=V*(i+a*j+stride*mode)+v;double value=mod(p);
     if(j){if(diag)value-=lower_((g*b+j)*a+i)*y(p-V*a);
      else for(int q=0;q<a;q++)value-=lower_(off+q*a+i)*y(V*(q+a*(j-1)+stride*mode)+v);}
     mod(p)=value;
    });team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,a),[&](int i){double value=0;for(int q=0;q<a;q++)value+=inv(off+q*a+i)*mod(V*(q+a*j+stride*mode)+v);y(V*(i+a*j+stride*mode)+v)=value;});team.team_barrier();
   }
   for(int j=b-1;j>=0;j--){int off=(g*b+j)*a*a;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,a),[&](int i){int p=V*(i+a*j+stride*mode)+v;double value=y(p);if(j+1<b)for(int q=0;q<a;q++)value-=T(off+q*a+i)*mod(V*(q+a*(j+1)+stride*mode)+v);mod(p)=value;});team.team_barrier();
   }
  });
 }else
#endif
 {
  Kokkos::parallel_for("modal inverse CPU",Range(0,N*V),KOKKOS_LAMBDA(int task){
   int mode=task/V,v=task%V,g=rows(task);
   bool borrowed=native.data()!=nullptr;
   const double*L_=borrowed?native(g).lu:L.data()+g*b*a*a;
   const double*T_=borrowed?native(g).transfer:T.data()+g*b*a*a;
   const double*lower_ptr=borrowed?native(g).lower:lower_.data()+g*b*(diag?a:a*a);
   const size_t*pivots=borrowed?native(g).permutation:nullptr;
   for(int j=0;j<b;j++){
    int off=j*a*a,poff=(g*b+j)*a;
    // First form every RHS before applying the pivot permutation.
    for(int i=0;i<a;i++){int p=V*(i+a*j+stride*mode)+v;double value=mod(p);
     if(j){if(diag)value-=lower_ptr[j*a+i]*y(p-V*a);
      else for(int q=0;q<a;q++)value-=lower_ptr[off+i*a+q]*y(V*(q+a*(j-1)+stride*mode)+v);}
     y(p)=value;
    }
    for(int i=0;i<a;i++){int pivot=borrowed?int(pivots[j*a+i]):perm(poff+i);mod(V*(i+a*j+stride*mode)+v)=y(V*(pivot+a*j+stride*mode)+v);}
    for(int i=0;i<a;i++){int p=V*(i+a*j+stride*mode)+v;double value=mod(p);for(int q=0;q<i;q++)value-=L_[off+i*a+q]*y(V*(q+a*j+stride*mode)+v);y(p)=value;}
    for(int i=a-1;i>=0;i--){int p=V*(i+a*j+stride*mode)+v;double value=y(p);for(int q=i+1;q<a;q++)value-=L_[off+i*a+q]*y(V*(q+a*j+stride*mode)+v);y(p)=value/L_[off+i*a+i];}
   }
   for(int j=b-1;j>=0;j--)for(int i=0;i<a;i++){int p=V*(i+a*j+stride*mode)+v,off=j*a*a;double value=y(p);
    if(j+1<b)for(int q=0;q<a;q++)value-=T_[off+i*a+q]*mod(V*(q+a*(j+1)+stride*mode)+v);mod(p)=value;
   }
  });
 }
 if(output_modal)copy(out,mod);
 else Kokkos::parallel_for("modal Fourier reconstruction",Range(0,total),KOKKOS_LAMBDA(int p){int point=p/V,v=p%V,k=point/stride,row=point%stride;double value=0;for(int mode=0;mode<N;mode++)value+=f(mode*N+k)*mod(V*(row+stride*mode)+v);out(p)=value;});
}
Modal import_BY_modal(TP_Modal*m){
 int shape[3];if(TP_modal_shape(m,shape))throw std::runtime_error("invalid BY modal shape");
 int na=shape[0],nb=shape[1],np=shape[2],total=na*nb*np,groups=np/2+1;
  const double*lu,*transfer,*lower,*scale,*fourier;const size_t*perm;
  if(TP_modal_factors(m,&lu,&transfer,&lower,&scale,&fourier,&perm))throw std::runtime_error("invalid BY modal factors");
  Modal M(na,nb,np,1,groups,false,false);size_t blocks=size_t(groups)*nb*na*na;
  if constexpr(Kokkos::SpaceAccessibility<Exec,Kokkos::HostSpace>::accessible){
   std::vector<ModalPointers>native(groups);size_t group_size=size_t(nb)*na*na,rows=size_t(nb)*na;
   for(int g=0;g<groups;g++)native[g]={lu+g*group_size,transfer+g*group_size,lower+g*group_size,perm+g*rows};
   M.native_blocks=decltype(M.native_blocks)("borrowed BY modal factors",groups);using H=Kokkos::View<const ModalPointers*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;Kokkos::deep_copy(M.native_blocks,H(native.data(),native.size()));
   M.forward=View(const_cast<double*>(fourier),np*np);
  }else{
   M.lu=upload(lu,blocks,"BY LU");M.transfer=upload(transfer,blocks,"BY transfer");M.lower=upload(lower,blocks,"BY lower");M.forward=upload(fourier,np*np,"BY Fourier");
  }
  std::vector<double>scaling(total);std::vector<int>permutation(size_t(groups)*nb*na),rows(np);
  for(size_t i=0;i<permutation.size();i++)permutation[i]=perm[i];
  for(int k=0;k<np;k++){int g=k<=np/2?k:k-np/2;rows[k]=g;for(int row=0;row<na*nb;row++)scaling[row+na*nb*k]=scale[row+na*nb*g];}
  M.row_scale=upload(scaling,"BY scaling");M.permutation=Indices("BY pivots",permutation.size());M.block=Indices("BY mode groups",np);
  using H=Kokkos::View<const int*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  Kokkos::deep_copy(M.permutation,H(permutation.data(),permutation.size()));Kokkos::deep_copy(M.block,H(rows.data(),rows.size()));M.prepare();
 return M;
}
}
/* Scoped to one Newton call; no owning Views outlive runtime finalization.
 * Fixed seed data is reused; masses are read again at every evaluation. */
struct Puncture_BYWorkspace {
 int n[3],total;puncture::Operator A;puncture::View geometry,input,rhs,solution,export_jet;
 std::vector<double> seed;
 static std::vector<double> current_seed(){
  std::vector<double> values{params_get_real("par_b")};
  for(const char*name:{"par_P_plus1","par_P_plus2","par_P_plus3","par_P_minus1","par_P_minus2","par_P_minus3","par_S_plus1","par_S_plus2","par_S_plus3","par_S_minus1","par_S_minus2","par_S_minus3"})values.push_back(params_get_real((char*)name));
  return values;
 }
 bool current()const{return seed==current_seed();}
 Puncture_BYWorkspace(const int*shape):n{shape[0],shape[1],shape[2]},total(n[0]*n[1]*n[2]),A(n,1,false),
  input("BY iterate",total),rhs("BY RHS",total),solution("BY solution",total),export_jet("BY export jet",total),seed(current_seed()){
  using namespace puncture;initialize_spectral(A,n,seed[0]);
  std::vector<double> geo(size_t(total)*4);
  double b=seed[0];
  for(int k=0;k<n[2];k++)for(int j=0;j<n[1];j++)for(int i=0;i<n[0];i++){
   int p=i+n[0]*(j+n[1]*k);double al=Pih*(2*i+1)/n[0],be=Pih*(2*j+1)/n[1],a=.5*(1-std::cos(al)),B=-std::cos(be),X=2*std::atanh(a),R=Pih+2*std::atan(B),phi=2*Pi*k/n[2];
   double x=b*std::cosh(X)*std::cos(R),r=b*std::sinh(X)*std::sin(R),y=r*std::cos(phi),z=r*std::sin(phi),sa=std::sin(al),sb=std::sin(be);
   geo[4*p]=std::sqrt((x-b)*(x-b)+y*y+z*z);geo[4*p+1]=std::sqrt((x+b)*(x+b)+y*y+z*z);geo[4*p+2]=BY_KKofxyz(x,y,z);geo[4*p+3]=sa*sb*sa*sb*sa*sb;
  }
  geometry=upload(geo,"BY immutable seed");A.coefficients=View("BY linear coefficients",size_t(total)*10);
 }
};

static int BY_capacity(int na,int nb,int np,const PK_Options*options){
 if(na<4||nb<4||np<4||na>256||nb>256||np>256||np%2||!options||options->restart<1||options->restart>4096)return PK_INVALID;
 const long double points=(long double)na*nb*np,blocks=(long double)(np/2+1)*nb*na*na;
 const long double krylov=options->method==PK_GMRES?2.L*options->restart+5:8;
 const long double small=16.L*(na*na+nb*nb+np*np)+64.L*na*nb+16.L*(np/2+1)*nb*na;
 // Includes fixed seed, export buffers, native retained jets/factors and
 // overlap of LU with the CUDA inverse. Conservative aggregate capacity.
 const long double device=8.L*points*(50+krylov)+32.L*blocks+small;
 const long double aggregate=device+24.L*blocks+8.L*points*80+small;
 int budget=params_get_int("TP_execution_memory_limit_mib");
 if(budget<16||budget>65536||aggregate>(long double)budget*1024*1024)return PK_ALLOCATION;
#ifdef KOKKOS_ENABLE_CUDA
 if constexpr(std::is_same_v<puncture::Exec,Kokkos::Cuda>){auto free=Puncture_execution_device_free_bytes();if(!free||device+512.L*1024*1024>free)return PK_ALLOCATION;}
#endif
 return PK_SUCCESS;
}
extern "C" {
const char *Puncture_execution_name(){return puncture::Exec::name();}
int Puncture_execution_initialize(int threads){
 puncture::ExecutionLock lock(puncture::execution_mutex());
 if(threads<0||Kokkos::is_finalized())return -1;
 auto&r=puncture::runtime();
 if(Kokkos::is_initialized())return threads&&Kokkos::DefaultHostExecutionSpace().concurrency()!=threads?-1:0;
 try{Kokkos::InitializationSettings settings;if(threads)settings.set_num_threads(threads);Kokkos::initialize(settings);r.owned=true;r.threads=Kokkos::DefaultHostExecutionSpace().concurrency();
  r.prior_callbacks=Kokkos::Tools::Experimental::get_callbacks();
  Kokkos::Tools::Experimental::set_allocate_data_callback(puncture::allocation);Kokkos::Tools::Experimental::set_deallocate_data_callback(puncture::deallocation);r.stats.memory_tracking_available=1;
  // Register AFTER Kokkos initializes its own static state, so teardown runs
  // before that state is destroyed. A Runtime destructor has the wrong order.
  std::atexit([](){auto&r=puncture::runtime();if(r.owned&&Kokkos::is_initialized()&&!Kokkos::is_finalized()){
   Kokkos::Tools::Experimental::set_allocate_data_callback(r.prior_callbacks.allocate_data);Kokkos::Tools::Experimental::set_deallocate_data_callback(r.prior_callbacks.deallocate_data);Kokkos::finalize();}});
  return threads&&r.threads!=threads?-1:0;}catch(...){return -1;}
}
int Puncture_execution_concurrency(){return Kokkos::is_initialized()?Kokkos::DefaultHostExecutionSpace().concurrency():0;}
unsigned long long Puncture_execution_device_free_bytes(){
#ifdef KOKKOS_ENABLE_CUDA
 if constexpr(std::is_same_v<puncture::Exec,Kokkos::Cuda>){size_t free=0,total=0;if(!Kokkos::is_initialized()||cudaMemGetInfo(&free,&total)!=cudaSuccess)return 0;return free;}
#endif
 return 0;
}
int Puncture_execution_device_description(char*buffer,int capacity){
 if(!buffer||capacity<1)return -1;
#ifdef KOKKOS_ENABLE_CUDA
 if constexpr(std::is_same_v<puncture::Exec,Kokkos::Cuda>){
  puncture::ExecutionLock lock(puncture::execution_mutex());int device=0,count=0;cudaDeviceProp prop{};size_t free=0,total=0;
  if(!Kokkos::is_initialized()||cudaGetDeviceCount(&count)!=cudaSuccess||cudaGetDevice(&device)!=cudaSuccess||cudaGetDeviceProperties(&prop,device)!=cudaSuccess||cudaMemGetInfo(&free,&total)!=cudaSuccess)return -1;
  char uuid[37];const auto*u=reinterpret_cast<const unsigned char*>(prop.uuid.bytes);
  snprintf(uuid,sizeof(uuid),"%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x",u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],u[9],u[10],u[11],u[12],u[13],u[14],u[15]);
  int used=snprintf(buffer,capacity,"{\"name\":\"%s\",\"uuid\":\"GPU-%s\",\"visible_ordinal\":%d,\"visible_count\":%d,\"total_bytes\":%llu,\"free_bytes\":%llu,\"sm_count\":%d,\"compute_major\":%d,\"compute_minor\":%d}",prop.name,uuid,device,count,(unsigned long long)total,(unsigned long long)free,prop.multiProcessorCount,prop.major,prop.minor);
  return used>=0&&used<capacity?0:-1;
 }
#endif
 return -1;
}
int Puncture_execution_statistics(PunctureExecutionStats*s){if(!s)return -1;puncture::ExecutionLock lock(puncture::execution_mutex());std::lock_guard<std::mutex>memory(puncture::runtime().memory_mutex);*s=puncture::statistics();return 0;}
void Puncture_execution_reset_statistics(){puncture::ExecutionLock lock(puncture::execution_mutex());std::lock_guard<std::mutex>memory(puncture::runtime().memory_mutex);auto&s=puncture::statistics();auto h=s.kokkos_host_current_bytes,d=s.kokkos_device_current_bytes;int enabled=s.memory_tracking_available;s={};s.kokkos_host_current_bytes=s.kokkos_host_peak_bytes=h;s.kokkos_device_current_bytes=s.kokkos_device_peak_bytes=d;s.memory_tracking_available=enabled;}
Puncture_BYWorkspace*Puncture_kokkos_BY_create(int na,int nb,int np){
 puncture::ExecutionLock lock(puncture::execution_mutex());
 try{
  puncture::initialize();PK_Options options{};options.method=params_get_int("TP_krylov_solver");options.restart=params_get_int("TP_krylov_restart");
  if(BY_capacity(na,nb,np,&options))return nullptr;
  auto seed=Puncture_BYWorkspace::current_seed();if(seed[0]<=0||!std::all_of(seed.begin(),seed.end(),[](double x){return std::isfinite(x);}))return nullptr;
  puncture::Timed timer(puncture::statistics().setup_seconds);const int shape[3]={na,nb,np};auto*w=new Puncture_BYWorkspace(shape);puncture::statistics().workspace_creations++;return w;
 }catch(const std::exception&e){fprintf(stderr,"Kokkos BY setup: %s\n",e.what());return nullptr;}
}
void Puncture_kokkos_BY_destroy(Puncture_BYWorkspace*workspace){puncture::ExecutionLock lock(puncture::execution_mutex());if(workspace){puncture::Exec().fence();delete workspace;}}
int Puncture_kokkos_BY_residual(Puncture_BYWorkspace*w,derivs*v,double*F,derivs*u){
 if(!w||!v||!F||!u||v==u||v->size!=w->total||u->size!=w->total||params_get_int("do_residuum_debug_output"))return PK_INVALID;
 const double*components[20]={v->d0,v->d1,v->d2,v->d3,v->d11,v->d12,v->d13,v->d22,v->d23,v->d33,u->d0,u->d1,u->d2,u->d3,u->d11,u->d12,u->d13,u->d22,u->d23,u->d33};
 if(!std::all_of(std::begin(components),std::end(components),[](const double*p){return p!=nullptr;}))return PK_INVALID;
 uintptr_t starts[21];starts[20]=reinterpret_cast<uintptr_t>(F);const uintptr_t bytes=sizeof(double)*w->total;
 for(int i=0;i<20;i++)starts[i]=reinterpret_cast<uintptr_t>(components[i]);
 for(int i=0;i<21;i++){if(starts[i]>UINTPTR_MAX-bytes)return PK_INVALID;for(int j=0;j<i;j++)if(starts[i]<starts[j]+bytes&&starts[j]<starts[i]+bytes)return PK_INVALID;}
 puncture::ExecutionLock lock(puncture::execution_mutex());
 try{
  using namespace puncture;if(!w->current())return PK_INVALID;Timed timer(statistics().nonlinear_seconds);statistics().nonlinear_calls++;
  {Timed timer(statistics().transfer_seconds);Kokkos::deep_copy(w->input,Host(v->d0,w->total));statistics().host_to_device_bytes+=8ull*w->total;}
  w->A.physical_fields(w->input);auto fields=w->A.fields,geo=w->geometry,out=w->rhs;double mp=params_get_real("par_m_plus"),mm=params_get_real("par_m_minus");
  Kokkos::parallel_for("BY nonlinear residual",Range(0,w->total),KOKKOS_LAMBDA(int p){
   double psi=1.+.5*mp/geo(4*p)+.5*mm/geo(4*p+1)+fields(10*p),psi2=psi*psi,psi4=psi2*psi2,psi7=psi*psi2*psi4;
   out(p)=(fields(10*p+4)+fields(10*p+7)+fields(10*p+9)+.125*geo(4*p+2)/psi7)*geo(4*p+3);
  });
  if(!finite(out)||!finite(fields))return PK_NONFINITE;download(out,F);
  double*raw[10]={v->d0,v->d1,v->d2,v->d3,v->d11,v->d12,v->d13,v->d22,v->d23,v->d33};
  double*cart[10]={u->d0,u->d1,u->d2,u->d3,u->d11,u->d12,u->d13,u->d22,u->d23,u->d33};
  for(int d=0;d<10;d++){
   if(d)download(w->A.spectral.work[d],raw[d]);
   auto jet=w->export_jet;Kokkos::parallel_for("BY export Cartesian jet",Range(0,w->total),KOKKOS_LAMBDA(int p){jet(p)=fields(10*p+d);});download(jet,cart[d]);
  }
  return PK_SUCCESS;
 }catch(const std::exception&e){fprintf(stderr,"Kokkos BY residual: %s\n",e.what());return PK_CALLBACK;}
}
int Puncture_kokkos_BY_linear(Puncture_BYWorkspace*w,const double*U,TP_Modal*m,
                       const double*rhs,double*solution,const PK_Options*options,PK_Result*result,PK_Monitor monitor,void*mc){
 if(result){*result={};result->status=PK_INVALID;result->true_residual=result->relative_residual=NAN;}
 if(!w||!U||!m||!rhs||!solution||!options||!result)return PK_INVALID;
 puncture::ExecutionLock lock(puncture::execution_mutex());
 int na=w->n[0],nb=w->n[1],np=w->n[2],shape[3];if(TP_modal_shape(m,shape)||shape[0]!=na||shape[1]!=nb||shape[2]!=np||!w->current())return result->status=PK_INVALID;
 try{
  int capacity=BY_capacity(na,nb,np,options);if(capacity)return result->status=capacity;const int total=w->total;
  auto&stats=puncture::statistics();double begin=puncture::seconds();
  auto&A=w->A;auto geo=w->geometry,coeff=A.coefficients,u=w->input;
  Kokkos::deep_copy(u,puncture::Host(U,total));stats.host_to_device_bytes+=8ull*total;
  double mp=params_get_real("par_m_plus"),mm=params_get_real("par_m_minus");
  Kokkos::parallel_for("BY freeze Jacobian",puncture::Range(0,total),KOKKOS_LAMBDA(int p){
   double psi=1.+.5*mp/geo(4*p)+.5*mm/geo(4*p+1)+u(p),p2=psi*psi,p4=p2*p2,weight=geo(4*p+3);
   for(int d=0;d<10;d++)coeff(10*p+d)=0;coeff(10*p)=-.875*geo(4*p+2)/(p4*p4)*weight;coeff(10*p+4)=coeff(10*p+7)=coeff(10*p+9)=weight;
  });
  auto M=puncture::import_BY_modal(m);
  auto v=w->rhs,x=w->solution;Kokkos::deep_copy(v,puncture::Host(rhs,total));Kokkos::deep_copy(x,puncture::Host(solution,total));stats.host_to_device_bytes+=16ull*total;
  puncture::Exec().fence();stats.setup_seconds+=puncture::seconds()-begin;
  stats.resident_bytes=std::max(stats.resident_bytes,8ull*(A.spectral.total*20+A.chain.extent(0)+A.coefficients.extent(0)+w->geometry.extent(0)+w->input.extent(0)+w->export_jet.extent(0)+M.lu.extent(0)+M.inverse.extent(0)+M.transfer.extent(0)+M.lower.extent(0)+M.workspace.extent(0)+M.column.extent(0)+M.row_scale.extent(0)+M.forward.extent(0)+v.extent(0)+x.extent(0)));
  // The native BY adapter supplies historical host dot/norm callbacks for
  // BiCGStab. Explicit Kokkos execution uses device L2 reductions instead.
  PK_Options device_options=*options;device_options.dot=nullptr;device_options.norm=nullptr;
  int status=puncture::solve(v,x,device_options,[&](auto in,auto out){A.apply(in,out);},[&](auto in,auto out){M.apply(in,out);},*result,monitor,mc);
  puncture::download(x,solution);return status;
 }catch(const std::exception&e){fprintf(stderr,"Kokkos BY: %s\n",e.what());return result->status=PK_CALLBACK;}
}
int Puncture_kokkos_BY(int na,int nb,int np,const double*U,TP_Modal*m,const double*rhs,double*solution,const PK_Options*options,PK_Result*result,PK_Monitor monitor,void*mc){
 if(result){*result={};result->status=PK_INVALID;result->true_residual=result->relative_residual=NAN;}
 if(!U||!m||!rhs||!solution||!options||!result)return PK_INVALID;
 auto*w=Puncture_kokkos_BY_create(na,nb,np);if(!w)return result->status=PK_ALLOCATION;
 int status=Puncture_kokkos_BY_linear(w,U,m,rhs,solution,options,result,monitor,mc);Puncture_kokkos_BY_destroy(w);return status;
}
}
