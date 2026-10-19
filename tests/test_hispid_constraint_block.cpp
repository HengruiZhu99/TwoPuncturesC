// Coupled, admissible manufactured correction on the actual production backend.
// This is a linear-algebra control, not a resolved physical wave experiment.
#include "../src/HiSpID_solver.cpp"
#include <cstdio>
int main(int argc,char**argv){
 const double separation=argc==2?std::stod(argv[1]):16.;
 if(argc>2||!std::isfinite(separation)||separation<=0)return 1;
 puncture::initialize();HiSpID_Config c;HiSpID_default_config(&c);
 c.n[0]=128;c.n[1]=96;c.n[2]=4;c.memory_limit_mib=8192;
 c.krylov_restart=32;c.max_krylov=240;c.far_radius=0;c.conformal_choice=1;
 c.inner_flatten=1;c.inner_min[0]=.1;c.inner_max[0]=.2;c.inner_min[1]=c.inner_max[1]=0;
 c.hole[0]={};c.hole[0].mass=1;c.hole[0].spin[2]=.9999;
 c.hole[1]={};c.hole[1].center[2]=-separation;
 std::unique_ptr<HiSpID_Data,HiSpID_Delete>s(HiSpID_create_with_seed_family(&c,HISPID_SEED_LIU,PUNCTURE_KOKKOS,0,0));
 if(!s){std::fprintf(stderr,"%s\n",HiSpID_last_error());return 1;}
 std::vector<double>known(s->ntotal,0.),rhs(s->ntotal),x,action(s->ntotal),residual(s->ntotal);
 for(int j=0;j<c.n[1];j++)for(int i=0;i<c.n[0];i++){
  const double sigma=.5*(1-std::cos(Pi*(i+.5)/c.n[0])),zeta=-std::cos(Pi*(j+.5)/c.n[1]);
  const double p=1e-8*sigma*sigma*std::pow(1-zeta*zeta,2);
  for(int v=0;v<4;v++)known[4*(i+c.n[0]*j)+v]=p*(v? .2*v:1.);
 }
 // Freeze the same zero base for all Jacobian applications.
 s->device->fields(s->values.data(),true);jvp(*s,known.data(),rhs.data());
 auto M=preconditioner(*s);double start=puncture::seconds();bool ok=linear_solve(*s,M,rhs,x,1e-7,PK_GMRES);
 double elapsed=puncture::seconds()-start;jvp(*s,x.data(),action.data());for(size_t i=0;i<rhs.size();i++)residual[i]=rhs[i]-action[i];
 double true_relative=norm2v(residual)/norm2v(rhs),error=0,scale=0;
 std::vector<Fields>a,b;fields(*s,known.data(),a,false);fields(*s,x.data(),b,false);
 for(size_t i=0;i<a.size();i++)for(int v=0;v<4;v++)for(int d=0;d<4;d++){
  error+=std::pow(a[i][v][d]-b[i][v][d],2);scale+=a[i][v][d]*a[i][v][d];
 }
 double field_error=std::sqrt(error/scale);
 std::printf("constraint-block coupled witness: passed=%d iterations=%d true_relative=%.17g field_relative=%.17g seconds=%.17g focus_separation=%.17g preconditioner=%s\n",ok,s->diag.krylov_iterations,true_relative,field_error,elapsed,separation,HiSpID_preconditioner_description(s.get()));
 return ok&&std::string(HiSpID_preconditioner_description(s.get()))=="constraint_block_upper_triangular"&&true_relative<=1e-7&&field_error<1e-4?0:1;
}
