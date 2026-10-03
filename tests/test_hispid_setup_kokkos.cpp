/* Tests the actual execution-space setup against retained host geometry and
 * GSL coordinates. No nonlinear solve or physical acceptance is implied. */
#include "../src/HiSpID_solver.cpp"
#include <cstdio>
#include <memory>
static int checks;
#define CHECK(c) do{checks++;if(!(c)){fprintf(stderr,"line%d: %s (%s)\n",__LINE__,#c,HiSpID_last_error());return 1;}}while(0)
template<class Real> HISPID_GEOMETRY_INLINE void pack_jet(const hispid::JetT<Real>&j,double*out,int&k){
 out[k++]=j.v;for(int a=0;a<4;a++)out[k++]=j.d[a];for(int a=0;a<4;a++)for(int b=0;b<4;b++)out[k++]=j.h[a][b];
}
template<class Real> HISPID_GEOMETRY_INLINE void pack_seed(const hispid::SeedT<Real>&s,double*out){
 int k=0;for(int a=0;a<3;a++)for(int b=0;b<3;b++){
  pack_jet(s.metric[a][b],out,k);pack_jet(s.physical[a][b],out,k);pack_jet(s.extrinsic[a][b],out,k);pack_jet(s.A[a][b],out,k);
 }pack_jet(s.psi,out,k);pack_jet(s.K,out,k);
}
HISPID_GEOMETRY_INLINE void pack_cache(const Cached&c,double*out){
 int k=0;out[k++]=c.psi;out[k++]=c.R;out[k++]=c.K;out[k++]=c.g;out[k++]=c.lapPsi;
 for(int a=0;a<3;a++){out[k++]=c.divM[a];out[k++]=c.gradK[a];}
 for(int a=0;a<9;a++){out[k++]=c.M[a];out[k++]=c.inv[a];}
 for(int d=0;d<10;d++){out[k++]=c.lap[d];out[k++]=c.far_correction[d];}
 for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<10;d++)out[k++]=c.vec[a][b][d];
 for(int a=0;a<9;a++)for(int b=0;b<3;b++)for(int d=0;d<4;d++)out[k++]=c.L[a][b][d];out[k++]=c.weight;
}
static double discrepancy(const std::vector<double>&a,const std::vector<double>&b,size_t*worst=nullptr){
 if(a.size()!=b.size())return INFINITY;double error=0;
 for(size_t i=0;i<a.size();i++){if(!std::isfinite(a[i])||!std::isfinite(b[i]))return INFINITY;double e=std::abs(a[i]-b[i])/(1+std::abs(a[i]));if(e>error){error=e;if(worst)*worst=i;}}return error;
}
/* diff() lowers a second-order jet to first order. Graph-slice K, A and
 * mean K therefore have valid values/gradients, but no complete Hessian.
 * Metric/psi remain second order; every exported slot must still be finite. */
static bool defined_seed_slot(size_t index){
 int jet=(index%(38*21))/21,slot=index%21;
 return slot<5 || (jet<36?jet%4<2:jet==36);
}
static int seed_controls(){
 constexpr int components=38*21;double worst=0,all_slot_worst=0;
 for(int kind=0;kind<6;kind++)for(int choice=0;choice<2;choice++){
  HiSpID_Hole h{};h.mass=1;
  if(kind==1){h.spin[0]=.57;h.spin[1]=.76;}
  if(kind==2){h.velocity[0]=.531;h.velocity[1]=.708;}
  if(kind==3){h.spin[0]=.2;h.spin[1]=.3;h.spin[2]=.4;h.velocity[0]=.3;h.velocity[1]=-.2;h.velocity[2]=.1;}
  if(kind==4)h.spin[2]=.99;
  if(kind==5)h.velocity[0]=std::sqrt(.99);
  double spin=hispid::point_norm<long double>(h.spin),throat=.5*std::sqrt(1-spin*spin),boost=kind==5?10:1;
  std::vector<double>points;const double directions[4][3]={{1,0,0},{0,1,0},{0,0,1},{.6,0,.8}};
  for(double r:{.03125,.125,throat,.75,3.,100.,1000.})for(const auto&dir:directions){points.push_back(r*dir[0]/boost);points.push_back(r*dir[1]);points.push_back(r*dir[2]);}
  const int count=points.size()/3;auto xyz=puncture::upload(points,"seed control points");puncture::View values("seed jet export",components*count);puncture::Indices status("seed control status",count);
  Kokkos::parallel_for("seed jet controls",puncture::Range(0,count),KOKKOS_LAMBDA(int p){hispid::SeedT<double>s;int code=hispid::seed_geometry(h,choice,xyz.data()+3*p,s,true);status(p)=code;if(!code)pack_seed(s,values.data()+components*p);});
  auto codes=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{},status);std::vector<double>actual(components*count),expected(actual.size());puncture::download(values,actual.data());
  for(int p=0;p<count;p++){CHECK(codes(p)==hispid::geometry_ok);hispid::Seed s;hispid::seed(h,choice,points.data()+3*p,s);pack_seed(s,expected.data()+components*p);}
  size_t index=0;double all_error=discrepancy(expected,actual,&index),error=0;all_slot_worst=std::max(all_slot_worst,all_error);
  CHECK(std::isfinite(all_error));size_t defined_index=0;
  for(size_t k=0;k<actual.size();k++)if(defined_seed_slot(k)){double e=std::abs(actual[k]-expected[k])/(1+std::abs(expected[k]));if(e>error){error=e;defined_index=k;}}
  worst=std::max(worst,error);fprintf(stderr,"seed kind%d choice%d defined derivatives %.3e all-slot diagnostic %.3e\n",kind,choice,error,all_error);
  if(all_error>1e-10){int point=index/components,jet=(index%components)/21,slot=index%21;
   const char*names[4]={"conformal_metric","physical_metric","extrinsic","Atilde"};
   fprintf(stderr,"worst diagnostic point%d xyz[%.17g,%.17g,%.17g] %s tensor[%d,%d] jet slot%d defined%d host%.17g execution%.17g\n",point,points[3*point],points[3*point+1],points[3*point+2],jet<36?names[jet%4]:jet==36?"psi":"meanK",jet/12,(jet/4)%3,slot,defined_seed_slot(index),expected[index],actual[index]);
  }
  if(error>1e-10){int point=defined_index/components,jet=(defined_index%components)/21,slot=defined_index%21;
   fprintf(stderr,"worst defined point%d xyz[%.17g,%.17g,%.17g] jet%d tensor[%d,%d] slot%d host%.17g execution%.17g\n",point,points[3*point],points[3*point+1],points[3*point+2],jet,jet/12,(jet/4)%3,slot,expected[defined_index],actual[defined_index]);
  }CHECK(error<=1e-10);
  if(kind>=4){
   // These are actual coefficients consumed by the equations, including
   // the Gamma10 throat where the unused partial-Hessian diagnostic failed.
   HiSpID_Config cfg;HiSpID_default_config(&cfg);cfg.hole[0]=h;cfg.hole[1].mass=0;
   cfg.conformal_choice=choice;cfg.omega[0]=cfg.omega[1]=0;
   constexpr int fields=248;puncture::View cached("extreme seed coefficient export",fields*count);
   Kokkos::parallel_for("extreme seed coefficient controls",puncture::Range(0,count),KOKKOS_LAMBDA(int p){
    hispid::BackgroundT<double>bg;int code=hispid::background_geometry(cfg,xyz.data()+3*p,bg,true);status(p)=code;
    if(!code){Cached c;hispid::cache(bg,c);pack_cache(c,cached.data()+fields*p);}
   });
   codes=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{},status);
   actual.resize(fields*count);expected.resize(actual.size());puncture::download(cached,actual.data());
   for(int p=0;p<count;p++){CHECK(codes(p)==hispid::geometry_ok);hispid::Background bg;
    hispid::background(cfg,points.data()+3*p,bg);Cached c;hispid::cache(bg,c);pack_cache(c,expected.data()+fields*p);
   }
   double cache_error=discrepancy(expected,actual,&index);
   fprintf(stderr,"extreme seed kind%d choice%d all coefficients %.3e\n",kind,choice,cache_error);
   if(cache_error>1e-10)fprintf(stderr,"worst coefficient point%zu field%zu host%.17g execution%.17g\n",index/fields,index%fields,expected[index],actual[index]);
   CHECK(cache_error<=1e-10);
  }
 }
 // Device status paths exclude singular points before any cache is built.
 puncture::Indices code("excluded puncture status",1);HiSpID_Hole h{};h.mass=1;
 Kokkos::parallel_for("excluded puncture",puncture::Range(0,1),KOKKOS_LAMBDA(int){double point[3]={};hispid::SeedT<double>s;code(0)=hispid::seed_geometry(h,0,point,s);});
 auto status=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{},code);CHECK(status(0)==hispid::geometry_puncture);
 printf("seed geometry defined-derivative maximum %.3e; all-slot diagnostic %.3e\n",worst,all_slot_worst);return 0;
}
int main(){
 if(std::numeric_limits<long double>::digits<=std::numeric_limits<double>::digits){
  fprintf(stderr,"SKIP: strict setup comparison requires a wider long-double reference; equal-precision legacy jets have cancellation errors. Use the x86 control build.\n");return 77;
 }
 puncture::initialize();CHECK(seed_controls()==0);
 using Context=std::unique_ptr<HiSpID_Data,decltype(&HiSpID_destroy)>;
 for(int choice=0;choice<2;choice++)for(int modified=0;modified<2;modified++){
  HiSpID_Config c;HiSpID_default_config(&c);c.n[0]=8;c.n[1]=12;c.n[2]=8;c.memory_limit_mib=2048;c.krylov_restart=16;
  c.conformal_choice=choice;c.inner_flatten=modified;c.far_radius=modified?8:0;
  c.hole[0].mass=.6;c.hole[1].mass=.4;c.hole[0].center[0]=3;c.hole[1].center[0]=-3;
  c.hole[0].spin[0]=.072;c.hole[0].spin[1]=.054;c.hole[0].spin[2]=.108;c.hole[1].spin[0]=-.032;c.hole[1].spin[1]=.048;c.hole[1].spin[2]=.016;
  c.hole[0].velocity[0]=.03;c.hole[0].velocity[1]=.06;c.hole[1].velocity[0]=-.02;c.hole[1].velocity[1]=-.07;
  c.inner_min[0]=c.inner_min[1]=.1;c.inner_max[0]=c.inner_max[1]=.3;
  if(choice&&modified){
   // Cyclic rotation and nonzero origin exercise lab-frame exports.
   for(auto&h:c.hole){
    double center[3]={h.center[2]+1,h.center[0]-2,h.center[1]+.5};
    double spin[3]={h.spin[2],h.spin[0],h.spin[1]},velocity[3]={h.velocity[2],h.velocity[0],h.velocity[1]};
    for(int a=0;a<3;a++){h.center[a]=center[a];h.spin[a]=spin[a];h.velocity[a]=velocity[a];}
   }
  }
  Context ref(HiSpID_create(&c),HiSpID_destroy),dev(HiSpID_create_with_geometry(&c,PUNCTURE_KOKKOS,1),HiSpID_destroy);
  CHECK(ref&&dev);CHECK(dev->geometry.empty());CHECK(dev->work==nullptr);CHECK(dev->device->geometry.extent(0)==size_t(dev->npt));
  for(int a=0;a<3;a++){CHECK(discrepancy(ref->derivatives.raw.D[a],dev->derivatives.raw.D[a])<=1e-12);CHECK(discrepancy(ref->derivatives.raw.D2[a],dev->derivatives.raw.D2[a])<=1e-12);}
  for(int a=0;a<2;a++){CHECK(discrepancy(ref->derivatives.coordinate[a],dev->derivatives.coordinate[a])<=1e-12);CHECK(discrepancy(ref->derivatives.first_map[a],dev->derivatives.first_map[a])<=1e-12);CHECK(discrepancy(ref->derivatives.second_map[a],dev->derivatives.second_map[a])<=1e-12);CHECK(discrepancy(ref->derivatives.coefficient[a],dev->derivatives.coefficient[a])<=1e-12);}
  CHECK(discrepancy(ref->derivatives.forward,dev->derivatives.forward)<=1e-12);CHECK(discrepancy(ref->derivatives.inverse_phi,dev->derivatives.inverse_phi)<=1e-12);CHECK(discrepancy(ref->derivatives.inverse_phi2,dev->derivatives.inverse_phi2)<=1e-12);
  puncture::Operator gsl(c.n,4,true);puncture::initialize_spectral(gsl,c.n,ref->b,ref->derivatives.coordinate.data());
  std::vector<double>a(gsl.chain.extent(0)),b(a.size());puncture::download(gsl.chain,a.data());puncture::download(dev->device->op.chain,b.data());CHECK(discrepancy(a,b)<=1e-11);
  constexpr int cache_components=248;static_assert(sizeof(Cached)==cache_components*sizeof(double));
  const int points=ref->npt;puncture::View packed("operator cache export",cache_components*points);auto geometry=dev->device->geometry;
  Kokkos::parallel_for("operator cache controls",puncture::Range(0,points),KOKKOS_LAMBDA(int p){pack_cache(geometry(p),packed.data()+cache_components*p);});
  a.resize(cache_components*points);b.resize(a.size());for(int p=0;p<points;p++)pack_cache(ref->geometry[p],a.data()+cache_components*p);puncture::download(packed,b.data());
  size_t worst_cache=0;double cache_error=discrepancy(a,b,&worst_cache);fprintf(stderr,"cache choice%d modified%d %.3e\n",choice,modified,cache_error);
  if(cache_error>1e-10){
   int p=worst_cache/cache_components,field=worst_cache%cache_components;double xyz[3];Fields empty{};
   transform(*ref,p%c.n[0],(p/c.n[0])%c.n[1],p/(c.n[0]*c.n[1]),empty,xyz);
   double comparisons[2];for(int stable=0;stable<2;stable++){hispid::BackgroundT<double>bg;Cached cached;double packed[cache_components];
    CHECK(hispid::background_geometry(ref->local,xyz,bg,stable)==hispid::geometry_ok);hispid::cache(bg,cached);pack_cache(cached,packed);comparisons[stable]=packed[field];}
   fprintf(stderr,"worst binary coefficient point%d xyz[%.17g,%.17g,%.17g] field%d host%.17g execution%.17g same-point-double-legacy%.17g stable%.17g\n",p,xyz[0],xyz[1],xyz[2],field,a[worst_cache],b[worst_cache],comparisons[0],comparisons[1]);
  }
  CHECK(cache_error<=1e-10);
  const int n=ref->ntotal;std::vector<double>v(n),direction(n),ar(n),br(n);for(int i=0;i<n;i++){v[i]=1e-7*std::cos(.23*i);direction[i]=1e-6*std::sin(.31*i);}
  CHECK(HiSpID_residual(ref.get(),v.data(),ar.data())==0);CHECK(HiSpID_residual(dev.get(),v.data(),br.data())==0);CHECK(discrepancy(ar,br)<=1e-10);
  CHECK(HiSpID_jvp(ref.get(),v.data(),direction.data(),ar.data())==0);CHECK(HiSpID_jvp(dev.get(),v.data(),direction.data(),br.data())==0);CHECK(discrepancy(ar,br)<=1e-10);
  auto reference_M=preconditioner(*ref),execution_M=preconditioner(*dev);CHECK(execution_M.col.empty());CHECK(reference_M.blocks.size()==execution_M.blocks.size());
  for(size_t block=0;block<reference_M.blocks.size();block++){CHECK(discrepancy(reference_M.blocks[block].lu,execution_M.blocks[block].lu)<=1e-9);CHECK(discrepancy(reference_M.blocks[block].transfer,execution_M.blocks[block].transfer)<=1e-9);}
  auto device_M=device_preconditioner(execution_M);auto rhs=puncture::upload(direction,"setup modal inverse RHS");puncture::View answer("setup modal inverse result",n);device_M.apply(rhs,answer);puncture::download(answer,br.data());reference_M.solve(direction.data(),ar.data());CHECK(discrepancy(ar,br)<=1e-9);
  CHECK(HiSpID_set_unknowns(ref.get(),v.data(),n)==0);CHECK(HiSpID_set_unknowns(dev.get(),v.data(),n)==0);make_coefficients(*ref);make_coefficients(*dev);CHECK(discrepancy(ref->coefficients,dev->coefficients)<=1e-12);
  std::vector<double>ra(9*points),rb(ra.size());
  CHECK(HiSpID_equation_samples(ref.get(),ra.data(),ra.data()+3*points,ra.data()+4*points,ra.data()+5*points)==0);
  CHECK(HiSpID_equation_samples(dev.get(),rb.data(),rb.data()+3*points,rb.data()+4*points,rb.data()+5*points)==0);CHECK(discrepancy(ra,rb)<=1e-10);CHECK(dev->basefields.empty());
  double sample_xyz[9]={.9,.6,-.4,4.1,.9,-.5,15,3,2},gradient_a[81],gradient_b[81];HiSpID_Point sample_a[3],sample_b[3];
  CHECK(HiSpID_sample_with_derivatives(ref.get(),3,sample_xyz,sample_a,gradient_a)==0);CHECK(HiSpID_sample_with_derivatives(dev.get(),3,sample_xyz,sample_b,gradient_b)==0);
  static_assert(sizeof(HiSpID_Point)==43*sizeof(double));a.resize(3*43);b.resize(a.size());std::memcpy(a.data(),sample_a,sizeof(sample_a));std::memcpy(b.data(),sample_b,sizeof(sample_b));CHECK(discrepancy(a,b)<=1e-10);
  a.assign(gradient_a,gradient_a+81);b.assign(gradient_b,gradient_b+81);CHECK(discrepancy(a,b)<=1e-10);
  HiSpID_SetupStatistics stats{int(sizeof(HiSpID_SetupStatistics))};CHECK(HiSpID_setup_statistics(dev.get(),&stats)==0);CHECK(stats.geometry_execution==1&&stats.scalar_digits==53);CHECK(stats.spectral_seconds>=0&&stats.geometry_seconds>=0&&stats.coefficient_seconds>=0);
  stats.struct_size=0;CHECK(HiSpID_setup_statistics(dev.get(),&stats)==-1);
 }
 HiSpID_Config c;HiSpID_default_config(&c);CHECK(!HiSpID_create_with_geometry(&c,PUNCTURE_REFERENCE,1));CHECK(!HiSpID_create_with_geometry(&c,PUNCTURE_KOKKOS,2));
 printf("HiSpID %s execution-space setup: %d checks\n",puncture::Exec::name(),checks);return 0;
}
