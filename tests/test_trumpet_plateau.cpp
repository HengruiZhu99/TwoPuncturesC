// The puncture representation relies on these exact flat plateau equations.
#include "HiSpID_internal.hpp"
#include <cmath>
#include <cstdio>
#include <algorithm>
int main() {
  HiSpID_Config c;HiSpID_default_config(&c);c.far_radius=0;c.inner_flatten=1;
  c.conformal_choice=0;
  double worst=0;int count=0;
  for(int h=0;h<2;h++) {
    auto &v=c.hole[h];v.mass=h?.4:.6;v.center[0]=h?-3:3;
    for(int d=0;d<3;d++){v.center[d]=d?0:v.center[d];v.spin[d]=v.mass*v.mass*(d==2?.3:.1);v.velocity[d]=(h?-1:1)*(.02+.01*d);}
    c.inner_min[h]=.1*v.mass;c.inner_max[h]=.3*v.mass;
  }
  hispid::Jet u[4];
  for(int k=0;k<4;k++) {
    u[k].v=.003*(k+1);
    for(int i=1;i<4;i++) {
      u[k].d[i]=.002*(k+i);
      for(int j=1;j<4;j++)u[k].h[i][j]=.001*(k+1)*(i+j);
    }
  }
  double expected[4]{};
  for(int j=1;j<4;j++)expected[0]+=u[0].h[j][j];
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)expected[i+1]+=u[i+1].h[j+1][j+1]+u[j+1].h[i+1][j+1]/3;
  for(int h=0;h<2;h++)for(int axis=0;axis<3;axis++) {
    double x[3];for(int d=0;d<3;d++)x[d]=c.hole[h].center[d]+(d==axis?.5*c.inner_min[h]:0);
    hispid::Background b;hispid::background(c,x,b,HISPID_SEED_TRUMPET_R0_M);
    if(b.g.v!=0)return 1;
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)if(b.opinv[i][j].v!=(i==j?1:0))return 2;
    double out[4];hispid::equations(b,u,out);
    for(int k=0;k<4;k++)worst=std::max(worst,std::abs(out[k]-expected[k]));
    count++;
  }
  std::printf("{\"points\":%d,\"flat_operator_absolute_error\":%.17g,\"passed\":%s}\n",count,worst,worst<1e-13?"true":"false");
  return worst<1e-13?0:3;
}
