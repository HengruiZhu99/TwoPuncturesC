#include "../src/PunctureModalProjection.hpp"
#include <iostream>
#ifdef PUNCTURES_KOKKOS
KOKKOS_INLINE_FUNCTION
#endif
double witness(int which){
 double x[4]={0,0,0,0},w[2]={1,-1},mean=0;
 if(which==0){double e=0x1p-27;x[0]=1-e;x[2]=1;w[0]=1+e;}
 if(which==1){x[0]=1;x[2]=1-0x1p-53;mean=-0x1p-53;}
 if(which==2){x[0]=x[2]=.375;mean=.375;}
 return puncture::compensated_projection(x,2,w,2,mean);
}
int main(int argc,char**argv){
 // Exact dyadic answers, independent of floating dot-product arithmetic:
 // (1+e)(1-e)-1=-e^2; shifted difference=2^-53; constant mode=0.
 const double expected[3]={-0x1p-54,0x1p-53,0};bool passed=true;
 for(int i=0;i<3;i++){double x=witness(i);passed&=x==expected[i];std::cout<<"host "<<i<<" "<<std::hexfloat<<x<<"\n";}
#ifdef PUNCTURES_KOKKOS
 Kokkos::initialize(argc,argv);
 {
  Kokkos::View<double*>result("projection witnesses",3);
  Kokkos::parallel_for("projection witnesses",3,KOKKOS_LAMBDA(int i){result(i)=witness(i);});
  auto host=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),result);
  for(int i=0;i<3;i++){passed&=host(i)==expected[i];std::cout<<"device "<<i<<" "<<std::hexfloat<<host(i)<<"\n";}
 }
 Kokkos::finalize();
#endif
 return passed?0:1;
}
