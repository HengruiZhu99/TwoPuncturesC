// Exact constant conformal-metric directions: cancellation must not create
// derivatives at either compactified end. Uses extended host jets.
#include "../src/HiSpID_geometry_kernels.hpp"
#include <cstdio>
#include <algorithm>
int main(){
 long double old_error=0,new_error=0;
 for(double velocity:{0.,.6})for(double radius:{1e-6,1e-3,.2,1.,1e3,1e8}){
  HiSpID_Hole h{};h.mass=1;h.velocity[0]=velocity;
  double x[3]={radius*.3,radius*.4,radius*.5};
  hispid::SeedT<long double>old,stable;
  if(hispid::seed_geometry(h,0,x,old,false)||hispid::seed_geometry(h,0,x,stable,true))return 1;
  // With boost along x, normalized yy and zz are exactly constant.
  for(int i:{1,2})for(int d=1;d<4;d++){
   old_error=std::max(old_error,std::abs(old.metric[i][i].d[d]));
   new_error=std::max(new_error,std::abs(stable.metric[i][i].d[d]));
   for(int e=1;e<4;e++){
    old_error=std::max(old_error,std::abs(old.metric[i][i].h[d][e]));
    new_error=std::max(new_error,std::abs(stable.metric[i][i].h[d][e]));
   }
  }
 }
 printf("constant transverse conformal metric: old derivative max %.17Lg, stable %.17Lg\n",old_error,new_error);
 return new_error==0?0:1;
}
