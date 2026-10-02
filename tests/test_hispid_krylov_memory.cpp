#include "../src/HiSpID_solver.cpp"
#include <iostream>
int main(){
 HiSpID_Config cfg;HiSpID_default_config(&cfg);cfg.n[0]=10;cfg.n[1]=14;cfg.n[2]=8;cfg.far_radius=0;
 cfg.hole[0].spin[2]=.03;cfg.hole[1].velocity[1]=-.025;
 HiSpID_Data*s=HiSpID_create(&cfg);if(!s)return 1;
 std::vector<double>res(s->ntotal),rhs(s->ntotal),a,b;residual(*s,s->values.data(),res.data());Sparse M=preconditioner(*s);
 int checks=0;
 for(int restart:{2,8,64})for(int limit:{1,9,35})for(bool zero:{false,true}){
  s->local.krylov_restart=restart;s->local.max_krylov=limit;
  for(int q=0;q<s->ntotal;q++)rhs[q]=zero?0:std::sin(.17*q)+std::cos(.29*q);
  s->diag.krylov_iterations=0;bool old=gmres(*s,M,rhs,a,1e-8,true);int count=s->diag.krylov_iterations;
  s->diag.krylov_iterations=0;bool now=gmres(*s,M,rhs,b,1e-8,false);
  if(old!=now||count!=s->diag.krylov_iterations||a.size()!=b.size()||std::memcmp(a.data(),b.data(),a.size()*sizeof(double))){std::cerr<<"GMRES mismatch\n";return 2;}checks++;
 }
 HiSpID_destroy(s);std::cout<<"eager/lazy GMRES bitwise controls: "<<checks<<" passed (zero RHS, early cap, restarts, exhausted limit)\n";return 0;
}
