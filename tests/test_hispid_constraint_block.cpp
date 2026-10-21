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
 HiSpID_Wave22 wave{int(sizeof(HiSpID_Wave22)),1,1,.001,.9760972454792889,4,12,0,2,4};
 if(HISPID_HALF_TURN_SOLVE){
  HiSpID_Data gate;gate.half_turn=true;
  if(HiSpID_set_axisymmetric(&gate,1)!=-1||gate.axisymmetric||HiSpID_set_axisymmetric(&gate,0)!=0)return 5;
  const double center[3]={0,0,0},point[3]={5,7,9},rotated[3]={-5,-7,9};const int sign[3]={-1,-1,1},ds[4]={1,-1,-1,1};
  hispid::Jet h[2][3][3],ht[2][3][3];hispid::wave22(wave,center,point,h[0],ht[0]);hispid::wave22(wave,center,rotated,h[1],ht[1]);
  long double error=0,norm=0;auto check=[&](long double a,long double b){error+=std::pow(a-b,2);norm+=a*a;};
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int time=0;time<2;time++){
   const auto&a=time?ht[0][i][j]:h[0][i][j];const auto&b=time?ht[1][i][j]:h[1][i][j];int t=sign[i]*sign[j];check(a.v,t*b.v);
   for(int d=0;d<4;d++){check(a.d[d],t*ds[d]*b.d[d]);for(int e=0;e<4;e++)check(a.h[d][e],t*ds[d]*ds[e]*b.h[d][e]);}
  }
  double relative=std::sqrt(double(error/norm));std::printf("TT metric and gdot jet half-turn covariance relative=%.17g\n",relative);std::fflush(stdout);if(relative>1e-13)return 3;
 }
 std::unique_ptr<HiSpID_Data,HiSpID_Delete>s(HISPID_HALF_TURN_SOLVE?HiSpID_create_with_wave22(&c,HISPID_SEED_LIU,PUNCTURE_KOKKOS,0,&wave):HiSpID_create_with_seed_family(&c,HISPID_SEED_LIU,PUNCTURE_KOKKOS,0,0));
 if(!s){std::fprintf(stderr,"%s\n",HiSpID_last_error());return 1;}
 std::vector<double>known(s->ntotal,0.),rhs(s->ntotal),x,action(s->ntotal),residual(s->ntotal);
 for(int j=0;j<c.n[1];j++)for(int i=0;i<c.n[0];i++){
  const double sigma=.5*(1-std::cos(Pi*(i+.5)/c.n[0])),zeta=-std::cos(Pi*(j+.5)/c.n[1]);
  const double p=1e-8*sigma*sigma*std::pow(1-zeta*zeta,2);
  if(HISPID_HALF_TURN_SOLVE){
   known[4*(i+c.n[0]*j)]=p;known[4*(i+c.n[0]*j)+1]=.2*p;
   known[4*(i+c.n[0]*(j+c.n[1]))+2]=.4*p;
   known[4*(i+c.n[0]*(j+3*c.n[1]))+3]=.6*p;
  }else for(int v=0;v<4;v++)known[4*(i+c.n[0]*j)+v]=p*(v? .2*v:1.);
 }
 if(HISPID_HALF_TURN_SOLVE){
  ::residual(*s,s->values.data(),residual.data());long double total=0,defect=0;const int N=c.n[0]*c.n[1],np=c.n[2];
  for(int k=0;k<np/2;k++)for(int q=0;q<4*N;q++){double a=residual[q+4*N*k],b=residual[q+4*N*(k+np/2)],e=a-(q%4>=2?-b:b);defect+=(long double)e*e;total+=(long double)a*a+(long double)b*b;}
  const double parity=total?std::sqrt(double(defect/total)):0.;std::printf("unprojected wave source half-turn parity relative=%.17g\n",parity);std::fflush(stdout);if(parity>1e-11)return 4;
 }
 // Freeze the same zero base for all Jacobian applications.
 s->device->fields(s->values.data(),true);jvp(*s,known.data(),rhs.data());
 if(HISPID_HALF_TURN_SOLVE){
  long double total=0,defect=0;const int N=c.n[0]*c.n[1],np=c.n[2];
  for(int k=0;k<np/2;k++)for(int q=0;q<4*N;q++){
   double a=rhs[q+4*N*k],b=rhs[q+4*N*(k+np/2)],e=a-(q%4>=2?-b:b);defect+=(long double)e*e;total+=(long double)a*a+(long double)b*b;
  }
  const double parity=std::sqrt(double(defect/total));std::printf("full Jacobian half-turn parity relative=%.17g\n",parity);std::fflush(stdout);
  if(parity>1e-11)return 2;
 }
 auto M=preconditioner(*s);double start=puncture::seconds();bool ok=linear_solve(*s,M,rhs,x,1e-7,PK_GMRES);
 double elapsed=puncture::seconds()-start;jvp(*s,x.data(),action.data());for(size_t i=0;i<rhs.size();i++)residual[i]=rhs[i]-action[i];
 double true_relative=norm2v(residual)/norm2v(rhs),error=0,scale=0;
 std::vector<Fields>a,b;fields(*s,known.data(),a,false);fields(*s,x.data(),b,false);
 for(size_t i=0;i<a.size();i++)for(int v=0;v<4;v++)for(int d=0;d<4;d++){
  error+=std::pow(a[i][v][d]-b[i][v][d],2);scale+=a[i][v][d]*a[i][v][d];
 }
 double field_error=std::sqrt(error/scale);
 std::printf("constraint-block coupled witness: passed=%d iterations=%d true_relative=%.17g field_relative=%.17g seconds=%.17g focus_separation=%.17g preconditioner=%s\n",ok,s->diag.krylov_iterations,true_relative,field_error,elapsed,separation,HiSpID_preconditioner_description(s.get()));
 return ok&&std::string(HiSpID_preconditioner_description(s.get()))==(HISPID_HALF_TURN_SOLVE?"constraint_block_upper_triangular_half_turn":"constraint_block_upper_triangular")&&true_relative<=1e-7&&field_error<1e-4?0:1;
}
