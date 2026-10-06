#include "HiSpID_trumpet.hpp"
// Standalone development probe: no elliptic context or public ABI changes.
extern "C" int trumpet_probe(const HiSpID_Hole *hole,int count,const double *xyz,
                              double *gamma,double *K,double *dg,double *dK) {
  for(int p=0;p<count;p++) {
    hispid::Seed s;
    int code=hispid::trumpet_seed_geometry(*hole,1,xyz+3*p,s);
    if(code)return code;
    for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
      gamma[9*p+3*i+j]=s.physical[i][j].v;
      K[9*p+3*i+j]=s.extrinsic[i][j].v;
      for(int d=0;d<3;d++) {
        dg[27*p+9*d+3*i+j]=s.physical[i][j].d[d+1];
        dK[27*p+9*d+3*i+j]=s.extrinsic[i][j].d[d+1];
      }
    }
  }
  return 0;
}
