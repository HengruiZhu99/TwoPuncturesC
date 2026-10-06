#include "HiSpID_trumpet.hpp"
#include "reference/kerr_trumpet_81c7d803.hpp"
#include <cstdio>
#include <algorithm>

int main() {
  double worst=0, minmargin=1e100, worstH=0, maxdt=0;
  int samples=0, failed=0;
  const double points[][3]={{.03125,0,0},{0,0,.2},{.3,-.2,.4},
                            {1,0,0},{3,4,-2},{60,-20,30}};
  for(double spin: {0.,.9,.99}) for(double speed: {0.,.5,std::sqrt(.99)}) {
    HiSpID_Hole hole{};hole.mass=1;hole.spin[2]=spin;hole.velocity[0]=speed;
    for(const auto &x:points) {
      hispid::Seed s;long double margin=0;
      int code=hispid::trumpet_seed_geometry(hole,1,x,s,&margin);
      ++samples; minmargin=std::min(minmargin,double(margin));
      if(code) {++failed;std::printf("rejected spin=%g speed=%g x=%g,%g,%g code=%d margin=%.17g\n",spin,speed,x[0],x[1],x[2],code,double(margin));continue;}
      if(speed==0) {
        kerr_trumpet::Geometry<double> ref;
        if(kerr_trumpet::Evaluate(1.,spin,x,ref)!=kerr_trumpet::Status::success)return 2;
        for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
          double e=std::abs(double(s.physical[i][j].v)-ref.gamma[i][j].value)/(1+std::abs(ref.gamma[i][j].value));
          worst=std::max(worst,e);
          e=std::abs(double(s.extrinsic[i][j].v)-ref.K[i][j].value)/(1+std::abs(ref.K[i][j].value));
          worst=std::max(worst,e);
        }
      }
      hispid::Jet inv[3][3],C[3][3][3];hispid::invert(s.physical,inv);
      hispid::connection(s.physical,inv,C);
      long double kk=0;
      for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
        maxdt=std::max(maxdt,std::abs(double(s.physical[i][j].d[0])));
        for(int k=0;k<3;k++)for(int l=0;l<3;l++)
          kk+=inv[i][k].v*inv[j][l].v*s.extrinsic[i][j].v*s.extrinsic[k][l].v;
      }
      const double H=double(hispid::curvature(inv,C)+s.K.v*s.K.v-kk);
      worstH=std::max(worstH,std::abs(H));
    }
  }
  std::printf("{\"samples\":%d,\"rejected\":%d,\"min_slice_margin\":%.17g,\"max_reference_scaled_difference\":%.17g,\"max_native_H\":%.17g,\"max_abs_metric_dt\":%.17g}\n",samples,failed,minmargin,worst,worstH,maxdt);
  // Native H is only an algebra smoke test; independent physical observer follows.
  return failed || worst>1e-11 || worstH>1e-7 || !(maxdt>0);
}
