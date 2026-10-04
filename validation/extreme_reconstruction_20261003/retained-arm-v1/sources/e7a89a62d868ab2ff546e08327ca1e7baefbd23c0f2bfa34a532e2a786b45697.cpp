// Fixed-state diagnostic only: private solver access, no solve or new public API.
#include "../src/HiSpID_solver.cpp"
#include "../tests/hispid_cartesian_oracle.hpp"
#include "hispid_checkpoint.hpp"  // Existing AthenaK portable reader, bound at build time.
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace {
void array(std::ostream&out,const double*x,int n){
 out<<'[';for(int i=0;i<n;i++){if(i)out<<',';if(!std::isfinite(x[i]))throw std::runtime_error("nonfinite diagnostic");out<<x[i];}out<<']';
}
void jets(std::ostream&out,const Fields&f,int count=10){
 out<<'[';for(int v=0;v<4;v++){if(v)out<<',';array(out,f[v].data(),count);}out<<']';
}
void physical(std::ostream&out,const HiSpID_Data&s,const double*x,const Fields&input){
 using hispid::Jet;
 hispid::Background bg;hispid::background(s.local,x,bg);Jet f[4];
 for(int v=0;v<4;v++){
  f[v].v=input[v][0];for(int d=1;d<4;d++)f[v].d[d]=input[v][d];
  for(int d=0;d<6;d++)f[v].h[hi[d]][hj[d]]=f[v].h[hj[d]][hi[d]]=input[v][4+d];
 }
 f[0]=f[0]+bg.far_correction;
 const Jet psi=bg.psi+f[0];
 if(!(psi.v>0))throw std::runtime_error("nonpositive continuous conformal factor");
 double equation[4];hispid::equations(bg,f,equation);
 equation[0]*=-8/std::pow(double(psi.v),5);
 for(int i=1;i<4;i++)equation[i]/=std::pow(double(psi.v),10);
 Jet L[3][3],metric[3][3],inv[3][3],K[3][3],C[3][3][3],trace;
 hispid::longitudinal(bg.metric,bg.C,f+1,L);
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  metric[i][j]=hispid::power(psi,4)*bg.metric[i][j];
  K[i][j]=(bg.M[i][j]+L[i][j])/(psi*psi)+metric[i][j]*bg.K/Jet(3);
 }
 hispid::invert(metric,inv);hispid::connection(metric,inv,C);
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)trace=trace+inv[i][j]*K[i][j];
 long double K2=0;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)for(int l=0;l<3;l++)
  K2+=inv[i][k].v*inv[j][l].v*K[i][j].v*K[k][l].v;
 const long double R=hispid::curvature(inv,C);
 Jet momentum[3][3];for(int i=0;i<3;i++)for(int j=0;j<3;j++)momentum[i][j]=K[i][j]-trace*metric[i][j];
 double HM[4]={double(R+trace.v*trace.v-K2),0,0,0};hispid::divergence(inv,C,momentum,HM+1);
 const double terms[3]={double(R),double(trace.v),double(K2)};
 out<<"{\"attenuation\":"<<double(bg.g.v)<<",\"psi\":"<<double(psi.v)<<",\"equation_physical\":";array(out,equation,4);
 out<<",\"tensor_constraints\":";array(out,HM,4);out<<",\"R_trace_K2\":";array(out,terms,3);out<<'}';
}
}

int main(int argc,char**argv){
 try{
  if(argc!=5)throw std::runtime_error("usage: reconstruct_iterate checkpoint source_sha points output_json");
  if(std::filesystem::exists(argv[4]))throw std::runtime_error("refuse to overwrite diagnostic output");
  const auto cp=hispid_import::Read(argv[1],argv[2]);
  if(cp.acceptance!="diagnostic")throw std::runtime_error("this diagnostic requires an explicitly failed checkpoint");
  hispid_import::Context owned(HiSpID_create_sampler(&cp.config));
  if(!owned)throw std::runtime_error(HiSpID_last_error());auto&s=*owned;
  if(HiSpID_set_unknowns(&s,cp.unknowns.data(),int(cp.unknowns.size())))throw std::runtime_error("unknown vector rejected");
  make_coefficients(s);
  // Only differentiation workspace is built: no background cache or Krylov state.
  s.derivatives.raw.initialize(s.local.n);
  std::array<std::vector<double>,10>storage;for(auto&v:storage)v.resize(s.ntotal);
  derivs w{};double**targets[10]={&w.d0,&w.d1,&w.d2,&w.d3,&w.d11,&w.d12,&w.d13,&w.d22,&w.d23,&w.d33};
  for(int d=0;d<10;d++)*targets[d]=storage[d].data();
  std::copy(cp.unknowns.begin(),cp.unknowns.end(),w.d0);s.derivatives.apply(4,&w);
  std::ifstream points(argv[3]);points.imbue(std::locale::classic());
  if(!points)throw std::runtime_error("cannot read witness points");
  int count;if(!(points>>count)||count<1||count>32)throw std::runtime_error("invalid point count");
  std::ostringstream out;out<<std::setprecision(17);
  out<<"{\"purpose\":\"fixed_iterate_continuous_Cartesian_polynomial_diagnostic\",\"solve_performed\":false,\"physical_acceptance\":false,\"parameterization\":\""<<HiSpID_unknown_parameterization()<<"\",\"residual_scaling\":\""<<HiSpID_residual_scaling()<<"\",\"scalar_digits\":"<<std::numeric_limits<long double>::digits<<",\"records\":[";
  for(int p=0;p<count;p++){
   int node;double lab[3],x[3];
   if(!(points>>node>>lab[0]>>lab[1]>>lab[2])||node<-1||node>=s.npt)throw std::runtime_error("invalid point row");
   for(double v:lab)if(!std::isfinite(v))throw std::runtime_error("nonfinite witness point");
   to_local(s,lab,x);Fields sampled{};sample_fields(s,x,sampled);
   const Fields AD=cartesian_polynomial(s,x);
   HiSpID_Point value{};double dg[27];
   if(HiSpID_sample_with_derivatives(&s,1,lab,&value,dg))throw std::runtime_error(HiSpID_last_error());
   if(p)out<<',';out<<"{\"source_node_index\":"<<node<<",\"lab_xyz\":";array(out,lab,3);
   out<<",\"local_xyz\":";array(out,x,3);out<<",\"Cartesian_jets\":";jets(out,AD);
   out<<",\"sampled_correction_value_gradient\":";jets(out,sampled,4);
   out<<",\"sampled_gamma\":";array(out,value.gamma,9);out<<",\"sampled_Kij\":";array(out,value.Kij,9);
   out<<",\"sampled_psi\":"<<value.psi<<",\"sampled_dgamma\":";array(out,dg,27);
   out<<",\"continuous\":";physical(out,s,x,AD);
   if(node>=0){
    const int i=node%s.local.n[0],j=(node/s.local.n[0])%s.local.n[1],k=node/(s.local.n[0]*s.local.n[1]);
    Fields grid{};double local[3];gather(&w,node,grid);transform(s,i,j,k,grid,local);
    out<<",\"collocation_local_xyz\":";array(out,local,3);
    out<<",\"collocation_jets\":";jets(out,grid);out<<",\"collocation\":";physical(out,s,local,grid);
   }
   out<<'}';std::cout<<"retained point "<<p<<" node "<<node<<" completed\n"<<std::flush;
  }
  std::string extra;if(points>>extra)throw std::runtime_error("extra witness point data");
  out<<"]}\n";std::ofstream result(argv[4]);result<<out.str();
  if(!result)throw std::runtime_error("cannot save diagnostic");
  return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
