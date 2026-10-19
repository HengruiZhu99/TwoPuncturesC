#include "HiSpID_internal.hpp"
#include "HiSpID_curved_preconditioner.hpp"
#include <cstdio>
int main(){
 HiSpID_Config c;HiSpID_default_config(&c);
 c.hole[0].mass=1;c.hole[0].center[0]=8;c.hole[0].center[1]=c.hole[0].center[2]=0;
 c.hole[0].spin[0]=.9999;c.hole[0].spin[1]=c.hole[0].spin[2]=0;
 c.hole[1].mass=0;c.hole[1].center[0]=-8;c.far_radius=0;c.omega[0]=c.omega[1]=0;c.inner_min[0]=.1;c.inner_max[0]=.2;
 for(double radius:{.12,.15,.17451721,.2,.25353545,.2939671,.35,1.}){
  double x[3]={8,radius,0};hispid::Background bg;hispid::background(c,x,bg,HISPID_SEED_LIU);hispid::Cached g;hispid::cache(bg,g);
  auto q=hispid::curved_symbol(g,8,x[0],x[1],x[2],0,0);
  printf("r %.9g mixed correlation %.9g\n",radius,q.mixed/(2*std::sqrt(q.second[0]*q.second[1])));
 }
}
