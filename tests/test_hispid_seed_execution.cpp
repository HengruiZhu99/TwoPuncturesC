/* Value-only public batch API, including exact-seed conventions, bounded
 * batches and device failures. No solve or binary acceptance is implied. */
#include "HiSpID.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>
static int checks;
#define CHECK(c) do{checks++;if(!(c)){fprintf(stderr,"line%d: %s (%s)\n",__LINE__,#c,HiSpID_last_error());return 1;}}while(0)
int main(){
 CHECK(Puncture_execution_initialize(0)==0);
 double worst=0;
 for(int kind=0;kind<6;kind++)for(int choice=0;choice<2;choice++){
  HiSpID_Hole h{};h.mass=1;h.center[0]=1;h.center[1]=-2;h.center[2]=.5;
  if(kind==1){h.spin[0]=.57;h.spin[1]=.76;}
  if(kind==2){h.velocity[0]=.531;h.velocity[1]=.708;}
  if(kind==3){h.spin[0]=.2;h.spin[1]=.3;h.spin[2]=.4;h.velocity[0]=.3;h.velocity[1]=-.2;h.velocity[2]=.1;}
  if(kind==4)h.spin[2]=.99;
  if(kind==5)h.velocity[0]=std::sqrt(.99);
  double xyz[18]={};const double displacement[6][3]={{.125,0,0},{0,.25,0},{0,0,.5},{.1,.2,.3},{1,-2,3},{100,200,-300}};
  for(int p=0;p<6;p++)for(int d=0;d<3;d++)xyz[3*p+d]=h.center[d]+displacement[p][d];
  HiSpID_Point ref[6],explicit_ref[6],device[6];
  CHECK(HiSpID_seed(&h,choice,6,xyz,ref)==0);
  CHECK(HiSpID_seed_with_execution(&h,choice,6,xyz,explicit_ref,PUNCTURE_REFERENCE)==0);
  CHECK(std::memcmp(ref,explicit_ref,sizeof(ref))==0);
  CHECK(HiSpID_seed_with_execution(&h,choice,6,xyz,device,PUNCTURE_KOKKOS)==0);
  static_assert(sizeof(HiSpID_Point)==43*sizeof(double));
  double a[6*43],b[6*43];std::memcpy(a,ref,sizeof(ref));std::memcpy(b,device,sizeof(device));
  double error=0;for(int k=0;k<6*43;k++){
   CHECK(std::isfinite(a[k])&&std::isfinite(b[k]));error=std::max(error,std::abs(a[k]-b[k])/(1+std::abs(a[k])));
  }
  worst=std::max(worst,error);CHECK(error<=1e-10);
  for(const auto&p:device){CHECK(p.attenuation==1);for(double correction:p.correction)CHECK(correction==0);}
 }
 // Cross a chunk boundary. Every point remains in lab coordinates, and a
 // reused caller output must be replaced by the fresh exact-seed values.
 HiSpID_Hole h{};h.mass=1;
 constexpr int n=4101;std::vector<double>xyz(3*n);std::vector<HiSpID_Point>out(n);
 for(int p=0;p<n;p++){xyz[3*p]=2;xyz[3*p+1]=3;xyz[3*p+2]=4;out[p].correction[0]=42;}
 HiSpID_Point one;CHECK(HiSpID_seed_with_execution(&h,0,1,xyz.data(),&one,PUNCTURE_KOKKOS)==0);
 PunctureExecutionStats before{},after{};CHECK(Puncture_execution_statistics(&before)==0);
 CHECK(HiSpID_seed_with_execution(&h,0,n,xyz.data(),out.data(),PUNCTURE_KOKKOS)==0);
 CHECK(Puncture_execution_statistics(&after)==0);
 CHECK(after.kokkos_host_current_bytes==before.kokkos_host_current_bytes);
 CHECK(after.kokkos_device_current_bytes==before.kokkos_device_current_bytes);
 for(const auto&p:out)CHECK(std::memcmp(&one,&p,sizeof(one))==0);
 double bad[6]={1,2,3,0,0,0};HiSpID_Point unused[2];
 CHECK(HiSpID_seed_with_execution(&h,0,2,bad,unused,PUNCTURE_KOKKOS)==-2);
 CHECK(std::strstr(HiSpID_last_error(),"point 1")!=nullptr);
 bad[0]=std::numeric_limits<double>::quiet_NaN();
 CHECK(HiSpID_seed_with_execution(&h,0,2,bad,unused,PUNCTURE_KOKKOS)==-1);
 CHECK(HiSpID_seed_with_execution(&h,0,2,bad,unused,2)==-1);
 h.velocity[0]=1;CHECK(HiSpID_seed_with_execution(&h,0,2,xyz.data(),unused,PUNCTURE_KOKKOS)==-1);
 printf("HiSpID %s public seed export: %d checks; value difference %.3e\n",Puncture_execution_name(),checks,worst);
 return 0;
}
