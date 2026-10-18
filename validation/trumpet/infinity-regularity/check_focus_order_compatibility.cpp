#include "HiSpID_focus.hpp"
#include "focus-c1-reference.hpp"
#include <fstream>
#include <cstdlib>
int main(int argc,char**argv){
 if(argc!=3||hispid::focus_rank!=2)return 1;int a=std::atoi(argv[1]),b=2*a;
 std::ofstream out(argv[2],std::ios::binary);double error=0;
 for(int m=0;m<=4;m++){
  auto old=hispid::reference_focus_weights(a,b,m),now=hispid::focus_weights(a,b,m);
  for(size_t p=0;p<old.size();p++)error=std::max(error,std::abs(old[p]-now[p]));
  out.write(reinterpret_cast<char*>(now.data()),now.size()*sizeof(double));
 }
 fprintf(stderr,"n=%d old/new row maximum difference %.17g\n",a,error);
 return out&&error<1e-15?0:1;
}
