#include "HiSpID_internal.hpp"
#include <algorithm>
#include <cstring>
#include <limits>
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
void background(const HiSpID_Config&c,const double*x,Background&b){geometry_error(background_geometry(c,x,b));}
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
 for(int k=0;k<3;k++)if(c.n[k]<4||c.n[k]>256||(k==2&&c.n[k]%2))return false;
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
 if(npt*bytes_per_point+block_bytes>(double)c.memory_limit_mib*1024*1024)return false;
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
  hispid::Seed s;hispid::seed(*h,choice,xyz+3*p,s);std::memset(out+p,0,sizeof(*out));
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   out[p].gamma[3*i+j]=s.physical[i][j].v;out[p].Kij[3*i+j]=s.extrinsic[i][j].v;
   out[p].conformal_metric[3*i+j]=s.metric[i][j].v;out[p].Atilde[3*i+j]=s.A[i][j].v;
   }
   out[p].psi=s.psi.v;out[p].mean_curvature=s.K.v;out[p].attenuation=1;
 }}catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
int HiSpID_operators(const HiSpID_Config*c,const double*x,const double*j,double*out){
 if(!c||!x||!j||!out||!hispid::valid(*c))return -1;
 try{
  hispid::Background b;hispid::background(*c,x,b);hispid::Jet u[4];
  const int ii[6]={1,1,1,2,2,3},jj[6]={1,2,3,2,3,3};
  for(int k=0;k<4;k++){u[k].v=j[10*k];for(int d=0;d<3;d++)u[k].d[d+1]=j[10*k+1+d];
   for(int d=0;d<6;d++)u[k].h[ii[d]][jj[d]]=u[k].h[jj[d]][ii[d]]=j[10*k+4+d];}
  out[0]=hispid::laplacian(b.inv,b.C,u[0]);hispid::Jet L[3][3];
  hispid::longitudinal(b.metric,b.C,u+1,L);hispid::divergence(b.inv,b.C,L,out+1);out[4]=b.R;
 }catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
}
