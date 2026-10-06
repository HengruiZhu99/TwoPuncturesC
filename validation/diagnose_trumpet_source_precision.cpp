// No solve: compare scalar source Fourier regularity before/after cache rounding.
#include "HiSpID_cache_kernels.hpp"
#include <cstdio>
#include <limits>
int main(){
 HiSpID_Config c{};c.conformal_choice=0;c.inner_flatten=1;c.attenuation_power=4;
 const double masses[2]={.6,.4},spins[2][3]={{.072,.054,.108},{-.032,.048,.016}},vel[2][3]={{.03,.06,.01},{-.02,-.07,.025}};
 for(int h=0;h<2;h++){c.hole[h].mass=masses[h];c.hole[h].center[0]=h?-3:3;c.omega[h]=.5;
  double s2=0,v2=0;for(int j=0;j<3;j++){c.hole[h].spin[j]=spins[h][j];c.hole[h].velocity[j]=vel[h][j];s2+=spins[h][j]*spins[h][j];v2+=vel[h][j]*vel[h][j];}
  double r=masses[h]*std::sqrt(1-s2/std::pow(masses[h],4))*std::sqrt(1-v2);c.inner_min[h]=.15*r;c.inner_max[h]=.5*r;
 }
 std::printf("# long_double_digits=%d; unweighted scalar source; Fourier amplitude 2/N\n",std::numeric_limits<long double>::digits);
 std::puts("x,rho,source_max,cached_m6,extended_m6,extended_rounded_m6,cached_extended_linf");
 const long double pi=std::acos(-1.L);constexpr int n=64;
 for(double x0:{0.,2.7,3.5})for(double rho:{.1,.03,.01,.003,.001,.0003,.0001}){
  long double cs[3]={},sn[3]={},maxsrc=0,maxerr=0;
  for(int k=0;k<n;k++){
   long double phi=2*pi*k/n;double x[3]={x0,rho*(double)std::cos(phi),rho*(double)std::sin(phi)};
   hispid::Background b{};int status=hispid::background_geometry(c,x,b,false,HISPID_SEED_TRUMPET_R0_M);if(status)return status;
   hispid::Cached q;hispid::cache(b,q);
   double a2=0;long double A2=0,lap=0;
   for(int i=0;i<3;i++)for(int j=0;j<3;j++){
    long double v=b.psi.h[i+1][j+1];for(int l=0;l<3;l++)v-=b.C[l][i][j].v*b.psi.d[l+1];lap+=b.inv[i][j].v*v;
    for(int l=0;l<3;l++)for(int m=0;m<3;m++){
     a2+=q.inv[3*i+l]*q.inv[3*j+m]*q.M[3*i+j]*q.M[3*l+m];
     A2+=b.inv[i][l].v*b.inv[j][m].v*b.M[i][j].v*b.M[l][m].v;
    }
   }
   double cached=q.g*(-q.psi*q.R/8-std::pow(q.psi,5)*q.K*q.K/12+a2/(8*std::pow(q.psi,7))+q.lapPsi);
   long double psi=b.psi.v,K=b.K.v,ext=b.g.v*(-psi*hispid::curvature(b.inv,b.C)/8-std::pow(psi,5)*K*K/12+A2/(8*std::pow(psi,7))+lap);
   long double v[3]={cached,ext,(double)ext};
   maxsrc=std::max(maxsrc,std::abs(ext));maxerr=std::max(maxerr,std::abs(v[0]-v[1]));
   for(int j=0;j<3;j++){cs[j]+=v[j]*std::cos(6*phi)*2/n;sn[j]+=v[j]*std::sin(6*phi)*2/n;}
  }
  std::printf("%.17g,%.17g,%.17Lg,%.17Lg,%.17Lg,%.17Lg,%.17Lg\n",x0,rho,maxsrc,std::hypot(cs[0],sn[0]),std::hypot(cs[1],sn[1]),std::hypot(cs[2],sn[2]),maxerr);
 }
}
