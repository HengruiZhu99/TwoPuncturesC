#include <algorithm>
#include <cmath>
#include <iostream>
#include "../src/HiSpID_cache_kernels.hpp"
int main(){
 HiSpID_Config c{};c.hole[0]={.6,{3,0,0},{.072,.054,.108},{.03,.06,.01}};
 c.hole[1]={.4,{-3,0,0},{-.032,.048,.016},{-.02,-.07,.025}};
 c.inner_min[0]=.0826630035142687;c.inner_min[1]=.0554766275110519;
 c.inner_max[0]=.468423686580856;c.inner_max[1]=.31436755589596077;
 c.omega[0]=c.omega[1]=.5;c.attenuation_power=4;c.inner_flatten=1;
 long double compatibility=0,old_defect=0;double exterior=0;int transitions=0;
 for(int family:{HISPID_SEED_QI,HISPID_SEED_TRUMPET_R0_M})
 for(double radius:{.04,.15,.27,.39,.6}){
  const double x[3]={3+.6*radius,.8*radius,0};hispid::Background b{};
  if(hispid::background_geometry(c,x,b,true,family))return 1;
  auto old=b;for(int k=0;k<3;k++)for(int i=0;i<3;i++)for(int j=0;j<3;j++)old.opC[k][i][j]=b.g*b.C[k][i][j];
  if(family==HISPID_SEED_QI){
   for(int k=0;k<3;k++)for(int i=0;i<3;i++)for(int j=0;j<3;j++)
    {
     const auto&a=old.opC[k][i][j];const auto&v=b.opC[k][i][j];
     if(a.v!=v.v)return 2;
     for(int d=0;d<4;d++){if(a.d[d]!=v.d[d])return 2;for(int e=0;e<4;e++)if(a.h[d][e]!=v.h[d][e])return 2;}
    }
   continue;
  }
  if(b.g.v>0&&b.g.v<1)transitions++;
  for(int d=0;d<3;d++)for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   long double defect=b.opmetric[i][j].d[d+1],previous=defect;
   for(int k=0;k<3;k++){
    defect-=b.opC[k][d][i].v*b.opmetric[k][j].v+b.opC[k][d][j].v*b.opmetric[i][k].v;
    previous-=old.opC[k][d][i].v*b.opmetric[k][j].v+old.opC[k][d][j].v*b.opmetric[i][k].v;
   }
   compatibility=std::max(compatibility,std::abs(defect)/(1+std::abs(b.opmetric[i][j].d[d+1])));
   old_defect=std::max(old_defect,std::abs(previous));
  }
  if(b.g.v==0||b.g.v==1){
   hispid::Cached a{},o{};hispid::cache(b,a);hispid::cache(old,o);
   for(int d=0;d<10;d++){
    exterior=std::max(exterior,std::abs(a.lap[d]-o.lap[d]));
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)exterior=std::max(exterior,std::abs(a.vec[i][j][d]-o.vec[i][j][d]));
   }
  }
 }
 std::cout<<"metric_compatibility="<<double(compatibility)<<" old_transition_defect="<<double(old_defect)
          <<" plateau_exterior_operator_difference="<<exterior<<" transition_points="<<transitions<<'\n';
 return compatibility<1e-15L&&old_defect>1e-3L&&exterior<1e-13&&transitions>=2?0:3;
}
