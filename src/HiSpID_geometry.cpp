#include "HiSpID_internal.hpp"
#include "HiSpID_tau.hpp"
#include <algorithm>
#include <cstring>
#include <limits>
#ifdef PUNCTURES_KOKKOS
#include "PunctureKokkos.hpp"
#endif
namespace hispid {
thread_local std::string last_error;
static long double norm(const double*a){return point_norm<long double>(a);}
static void geometry_error(int status){
 if(status==geometry_puncture)throw std::runtime_error("sampling at unsmoothed puncture");
 if(status==geometry_determinant)throw std::runtime_error("nonpositive spatial determinant");
 if(status==geometry_nonspacelike)throw std::runtime_error("boosted slice is not spacelike");
 if(status==geometry_nonfinite)throw std::runtime_error("nonfinite background geometry");
}
void seed(const HiSpID_Hole&h,int choice,const double*x,Seed&s){geometry_error(seed_geometry(h,choice,x,s));}
void background(const HiSpID_Config&c,const double*x,Background&b,int family){geometry_error(background_geometry(c,x,b,false,family));}
double laplacian(const Jet inv[3][3],const Jet C[3][3][3],const Jet&u){return geometry_laplacian(inv,C,u);}

void longitudinal(const Jet metric[3][3],const Jet C[3][3][3],const Jet b[3],Jet L[3][3]){
 Jet db[3][3],div=0;for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  db[i][j]=diff(b[j],i+1);for(int k=0;k<3;k++)db[i][j]=db[i][j]+C[j][i][k]*b[k];
  if(i==j)div=div+db[i][j];
 }
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  L[i][j]=-Jet(2.0L/3)*metric[i][j]*div;
  for(int k=0;k<3;k++)L[i][j]=L[i][j]+metric[i][k]*db[j][k]+metric[j][k]*db[i][k];
 }
}

void divergence(const Jet inv[3][3],const Jet C[3][3][3],const Jet T[3][3],double out[3]){
 Jet up[3][3];for(int i=0;i<3;i++)for(int j=0;j<3;j++)
  for(int k=0;k<3;k++)for(int l=0;l<3;l++)up[i][j]=up[i][j]+inv[i][k]*inv[j][l]*T[k][l];
 for(int i=0;i<3;i++){
  long double v=0;for(int j=0;j<3;j++){
   v+=up[i][j].d[j+1];for(int k=0;k<3;k++)v+=C[i][j][k].v*up[k][j].v+C[j][j][k].v*up[i][k].v;
  }out[i]=(double)v;
 }
}

static long double contract(const Background&b,const Jet A[3][3],const Jet D[3][3]){
 long double q=0;for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)for(int l=0;l<3;l++)
  q+=b.inv[i][k].v*b.inv[j][l].v*A[i][j].v*D[k][l].v;
 return q;
}

void equations(const Background&b,const Jet u[4],double out[4],const Jet *du){
 Jet L[3][3],A[3][3];longitudinal(b.metric,b.C,u+1,L);
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)A[i][j]=b.M[i][j]+L[i][j];
 long double psi=b.psi.v+u[0].v;
 if(!(psi>0)){for(int k=0;k<4;k++)out[k]=std::numeric_limits<double>::quiet_NaN();return;}
 long double A2=contract(b,A,A),K=b.K.v,g=b.g.v;
 if(!du){
  out[0]=laplacian(b.opinv,b.opC,u[0])+g*(-psi*b.R/8-std::pow(psi,5)*K*K/12+A2/(8*std::pow(psi,7))+b.lapPsi);
  Jet OL[3][3];longitudinal(b.opmetric,b.opC,u+1,OL);
  divergence(b.opinv,b.opC,OL,out+1);
  for(int i=0;i<3;i++){
   long double grad=0;for(int j=0;j<3;j++)grad+=b.inv[i][j].v*b.K.d[j+1];
   out[i+1]+=g*(b.divM[i]-2.0L/3*std::pow(psi,6)*grad);
  }
 }else{
  Jet DL[3][3],OL[3][3];longitudinal(b.metric,b.C,du+1,DL);
  out[0]=laplacian(b.opinv,b.opC,du[0])+g*((-b.R/8-5*std::pow(psi,4)*K*K/12-7*A2/(8*std::pow(psi,8)))*du[0].v
     +contract(b,A,DL)/(4*std::pow(psi,7)));
  longitudinal(b.opmetric,b.opC,du+1,OL);divergence(b.opinv,b.opC,OL,out+1);
  for(int i=0;i<3;i++){
   long double grad=0;for(int j=0;j<3;j++)grad+=b.inv[i][j].v*b.K.d[j+1];
   out[i+1]-=g*4*std::pow(psi,5)*du[0].v*grad;
  }
 }
}

bool valid(const HiSpID_Config&c,bool sampler_only){
 if(c.conformal_choice<0||c.conformal_choice>1||c.attenuation_power<2||c.attenuation_power%2)return false;
 if(!std::isfinite(c.far_radius)||!std::isfinite(c.tolerance)||c.tolerance<=0||c.max_newton<0||c.max_krylov<1||c.krylov_restart<2||c.krylov_restart>200)return false;
 const int limit[3]={HISPID_MAX_RADIAL_POINTS,HISPID_MAX_POLAR_POINTS,HISPID_MAX_AZIMUTHAL_POINTS};
 for(int k=0;k<3;k++)if(c.n[k]<4||c.n[k]>limit[k]||(k==2&&c.n[k]%2))return false;
 /* Bound the compact cache, four-field modal FD stencil and Krylov basis before
  * allocation. Larger grids require an explicit per-context budget. */
 if(c.memory_limit_mib<16||c.memory_limit_mib>65536)return false;
 const double npt=(double)c.n[0]*c.n[1]*c.n[2];
 // Sampler: values, coefficients and SpecCoef's two padded scalar scratch
 // arrays/pointer tables, with headroom for the smallest allowed grids.
 // One16-byte azimuthal average per meridional point. Newton vector reuse
 // moves the existing factor bank and does not allocate a second bank.
 const double bytes_per_point=sampler_only?128:3192+4*(5*12+48)+32*(2*c.krylov_restart+30)+16.0/c.n[2];
 // Two dense matrices per radial block, scalar/vector factor groups, shared
 // by real-Fourier partners. Add permutations/diagonal couplings explicitly.
 const double block_bytes=sampler_only?0:2.0*(c.n[2]/2+1)*c.n[1]*(16.0*c.n[0]*c.n[0]+24.0*c.n[0]);
 double extra=sampler_only?0:double(polar_border_bytes(c.n[0],c.n[1],c.n[2]));
#if HISPID_AXIS_TAU
 if(!sampler_only)extra+=8.0*((2*c.n[0]+c.n[1]-2)*4*c.n[2]+2*c.n[1]+c.n[0]);
#endif
#if HISPID_STABLE_SCALAR_SOURCE
 if(!sampler_only)extra+=16*npt;
#endif
 if(npt*bytes_per_point+block_bytes+extra>(double)c.memory_limit_mib*1024*1024)return false;
 bool active=false;for(int h=0;h<2;h++){
  const auto&v=c.hole[h];if(!std::isfinite(v.mass)||v.mass<0)return false;
  for(int k=0;k<3;k++)if(!std::isfinite(v.center[k])||!std::isfinite(v.spin[k])||!std::isfinite(v.velocity[k]))return false;
  if(v.mass>0){active=true;if(norm(v.spin)>=v.mass*v.mass||norm(v.velocity)>=1)return false;}
  if(!std::isfinite(c.omega[h])||!std::isfinite(c.inner_min[h])||!std::isfinite(c.inner_max[h]))return false;
  if(c.inner_max[h]>0&&(c.inner_min[h]<0||c.inner_min[h]>=c.inner_max[h]))return false;
 }return active;
}
}

extern "C" {
const char *HiSpID_last_error(void){return hispid::last_error.c_str();}
void HiSpID_default_config(HiSpID_Config*c){
 if(!c)return;std::memset(c,0,sizeof(*c));
 c->hole[0].mass=c->hole[1].mass=.5;c->hole[0].center[0]=3;c->hole[1].center[0]=-3;
 c->n[0]=c->n[1]=16;c->n[2]=8;c->conformal_choice=1;c->inner_flatten=1;
 c->omega[0]=c->omega[1]=.5;c->attenuation_power=4;
 for(int h=0;h<2;h++){c->inner_min[h]=.05;c->inner_max[h]=.1;}
 c->far_radius=40;c->tolerance=1e-10;c->max_newton=12;c->max_krylov=600;c->krylov_restart=40;
 c->memory_limit_mib=2048;
}
int HiSpID_seed(const HiSpID_Hole*h,int choice,int count,const double*xyz,HiSpID_Point*out){
 if(!h||!xyz||!out||count<0||!(h->mass>0)||choice<0||choice>1)return -1;
 HiSpID_Config c;HiSpID_default_config(&c);c.hole[0]=*h;c.hole[1].mass=0;
 if(!hispid::valid(c))return -1;
 try{for(int p=0;p<count;p++){
  hispid::Seed s;hispid::seed(*h,choice,xyz+3*p,s);hispid::seed_values(s,out[p]);
 }}catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
int HiSpID_seed_with_execution(const HiSpID_Hole*h,int choice,int count,const double*xyz,HiSpID_Point*out,int execution){
 return HiSpID_seed_with_family(h,choice,count,xyz,out,execution,HISPID_SEED_QI);
}
int HiSpID_seed_with_family(const HiSpID_Hole*h,int choice,int count,const double*xyz,HiSpID_Point*out,int execution,int family){
 if(family!=HISPID_SEED_QI && family!=HISPID_SEED_TRUMPET_R0_M){hispid::last_error="invalid seed family";return -1;}
 if(execution==PUNCTURE_REFERENCE && family==HISPID_SEED_QI)return HiSpID_seed(h,choice,count,xyz,out);
 if(execution==PUNCTURE_REFERENCE){
  if(!h||!xyz||!out||count<0||!(h->mass>0)||choice<0||choice>1)return -1;
  HiSpID_Config c;HiSpID_default_config(&c);c.hole[0]=*h;c.hole[1].mass=0;
  if(!hispid::valid(c))return -1;
  try{for(int p=0;p<count;p++){
   for(int d=0;d<3;d++)if(!std::isfinite(xyz[3*p+d]))throw std::runtime_error("nonfinite seed coordinate");
   hispid::Seed s;hispid::geometry_error(hispid::seed_geometry(*h,choice,xyz+3*p,s,false,family));
   hispid::seed_values(s,out[p]);
   if(!hispid::finite_seed_values(out[p]))throw std::runtime_error("nonfinite seed output");
  }}catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
 }
 if(execution!=PUNCTURE_KOKKOS||!h||!xyz||!out||count<0||!(h->mass>0)||choice<0||choice>1){
  hispid::last_error="invalid execution seed request";return -1;
 }
 HiSpID_Config cfg;HiSpID_default_config(&cfg);cfg.hole[0]=*h;cfg.hole[1].mass=0;
 if(!hispid::valid(cfg)){hispid::last_error="invalid execution seed parameters";return -1;}
 for(size_t k=0;k<size_t(count)*3;k++)if(!std::isfinite(xyz[k])){hispid::last_error="nonfinite execution seed coordinate";return -1;}
#ifdef PUNCTURES_KOKKOS
 try{
  puncture::ExecutionLock lock(puncture::execution_mutex());puncture::initialize();
  const HiSpID_Hole hole=*h;
  // Bound temporary export memory for independent Cartesian FD/quadrature
  // batches. No retained host or device jet bank is needed.
  for(int first=0;first<count;){
   const int n=std::min(count-first,4096);
   auto x=puncture::upload(xyz+size_t(3)*first,size_t(3)*n,"execution seed coordinates");
   Kokkos::View<HiSpID_Point*,puncture::Exec>values(Kokkos::view_alloc(Kokkos::WithoutInitializing,"execution seed values"),n);
   puncture::Indices codes("execution seed status",n);
   Kokkos::parallel_for("execution seed value export",puncture::Range(0,n),KOKKOS_LAMBDA(int p){
    hispid::SeedT<double>s;int code=hispid::seed_geometry(hole,choice,x.data()+3*p,s,true,family);
    if(!code){hispid::seed_values(s,values(p));if(!hispid::finite_seed_values(values(p)))code=hispid::geometry_nonfinite;}
    codes(p)=code;
   });
   int failed=n;
   Kokkos::parallel_reduce("first failed execution seed",puncture::Range(0,n),KOKKOS_LAMBDA(int p,int&index){if(codes(p)&&p<index)index=p;},Kokkos::Min<int>(failed));
   if(failed<n){int code;Kokkos::deep_copy(code,Kokkos::subview(codes,failed));
    throw std::runtime_error("execution seed failed at point "+std::to_string(first+failed)+" (status "+std::to_string(code)+")");}
   using HostPoints=Kokkos::View<HiSpID_Point*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
   {puncture::Timed timer(puncture::statistics().transfer_seconds);
    Kokkos::deep_copy(HostPoints(out+first,n),values);puncture::statistics().device_to_host_bytes+=size_t(n)*sizeof(HiSpID_Point);}
   first+=n;
  }
 }catch(const std::exception&e){hispid::last_error=e.what();return -2;}
 return 0;
#else
 hispid::last_error="library was built without Kokkos";return -1;
#endif
}
int HiSpID_operators(const HiSpID_Config*c,const double*x,const double*j,double*out){
 return HiSpID_operators_with_seed_family(c,x,j,out,HISPID_SEED_QI);
}
int HiSpID_operators_with_seed_family(const HiSpID_Config*c,const double*x,const double*j,double*out,int family){
 if(family!=HISPID_SEED_QI && family!=HISPID_SEED_TRUMPET_R0_M)return -1;
 if(!c||!x||!j||!out||!hispid::valid(*c))return -1;
 try{
  hispid::Background b;hispid::background(*c,x,b,family);hispid::Jet u[4];
  const int ii[6]={1,1,1,2,2,3},jj[6]={1,2,3,2,3,3};
  for(int k=0;k<4;k++){u[k].v=j[10*k];for(int d=0;d<3;d++)u[k].d[d+1]=j[10*k+1+d];
   for(int d=0;d<6;d++)u[k].h[ii[d]][jj[d]]=u[k].h[jj[d]][ii[d]]=j[10*k+4+d];}
  out[0]=hispid::laplacian(b.inv,b.C,u[0]);hispid::Jet L[3][3];
  hispid::longitudinal(b.metric,b.C,u+1,L);hispid::divergence(b.inv,b.C,L,out+1);out[4]=b.R;
 }catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
}
