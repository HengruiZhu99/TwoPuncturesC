/* Actual execution-space fields, coupled Jacobian and conditioned production
 * modal inverse, compared to the independent reference implementation. */
#include "../src/HiSpID_solver.cpp"
#include <cstdio>
#include <memory>
static int checks;
#define CHECK(c) do{checks++;if(!(c)){fprintf(stderr,"line%d: %s\n",__LINE__,#c);return 1;}}while(0)
static double relative(const std::vector<double>&a,const std::vector<double>&b){double d=0,scale=0;for(size_t i=0;i<a.size();i++){d+=(a[i]-b[i])*(a[i]-b[i]);scale+=a[i]*a[i];}return std::sqrt(d)/(1e-13+std::sqrt(scale));}
int main(){
 puncture::initialize();HiSpID_Config c;HiSpID_default_config(&c);c.n[0]=12;c.n[1]=24;c.n[2]=16;c.memory_limit_mib=2048;c.krylov_restart=16;
 c.hole[0].mass=.6;c.hole[1].mass=.4;c.hole[0].center[0]=3;c.hole[1].center[0]=-3;
 double spins[2][3]={{.072,.054,.108},{-.032,.048,.016}},vel[2][3]={{.03,.06,.01},{-.02,-.07,.025}};
 for(int h=0;h<2;h++)for(int d=0;d<3;d++){c.hole[h].spin[d]=spins[h][d];c.hole[h].velocity[d]=vel[h][d];}
 std::unique_ptr<HiSpID_Data,decltype(&HiSpID_destroy)>ref(HiSpID_create(&c),HiSpID_destroy),dev(HiSpID_create_with_execution(&c,1),HiSpID_destroy);
 CHECK(ref&&dev);const int n=ref->ntotal;std::vector<double>v(n),d(n),ar(n),br(n);
 for(int i=0;i<n;i++){v[i]=1e-7*std::cos(.23*i);d[i]=1e-6*std::sin(.31*i);}
 std::vector<Fields>af,bf;fields(*ref,v.data(),af,false);fields(*dev,v.data(),bf,false);
 std::vector<double>fa(af.size()*40),fb(fa.size());std::memcpy(fa.data(),af.data(),fa.size()*8);std::memcpy(fb.data(),bf.data(),fb.size()*8);
 double fields_error=relative(fa,fb);CHECK(fields_error<=1e-10);
 for(size_t i=0;i<fa.size();i++)CHECK(std::abs(fa[i]-fb[i])/(1+std::abs(fa[i]))<=1e-10);
 CHECK(HiSpID_residual(ref.get(),v.data(),ar.data())==0);CHECK(HiSpID_residual(dev.get(),v.data(),br.data())==0);CHECK(relative(ar,br)<=1e-10);
 CHECK(HiSpID_jvp(ref.get(),v.data(),d.data(),ar.data())==0);CHECK(HiSpID_jvp(dev.get(),v.data(),d.data(),br.data())==0);double jvp_error=relative(ar,br);CHECK(jvp_error<=1e-10);
 auto M=preconditioner(*ref);auto compact=preconditioner(*dev);CHECK(compact.col.empty()&&compact.val.empty());CHECK(M.blocks.size()==compact.blocks.size());
 for(size_t group=0;group<M.blocks.size();group++){CHECK(relative(M.blocks[group].lu,compact.blocks[group].lu)<=1e-9);CHECK(relative(M.blocks[group].transfer,compact.blocks[group].transfer)<=1e-9);}
 auto deviceM=device_preconditioner(compact);auto rhs=puncture::upload(d,"conditioned RHS");puncture::View x("inverse result",n);deviceM.apply(rhs,x);puncture::download(x,br.data());M.solve(d.data(),ar.data());double inverse_error=relative(ar,br);CHECK(inverse_error<=1e-9);
 // Use a known bounded modal unknown to manufacture a balanced physical
 // RHS. Arbitrary physical noise is amplified by tiny high-mode row scales;
 // that conditioning floor is retained above in the inverse comparison.
 std::vector<double>known(n),balanced(n),weighted(n,0);
 for(int row=0;row<n;row++)known[row]=std::sin(.19*row);
 for(int row=0;row<n;row++){for(size_t q=0;q<M.col[row].size();q++)weighted[row]+=M.val[row][q]*known[M.col[row][q]];weighted[row]*=M.row_scale[row];}
 ref->derivatives.raw.along(2,ref->derivatives.inverse,4,weighted.data(),balanced.data());
 auto balanced_rhs=puncture::upload(balanced,"manufactured modal RHS");deviceM.apply(balanced_rhs,x);puncture::download(x,br.data());
 // Independently apply the original weighted FD stencil to the returned
 // unknown, then reconstruct physical Fourier values; no inverse code used.
 std::vector<double>modal(n,0),physical(n);
 for(int row=0;row<n;row++){for(size_t q=0;q<M.col[row].size();q++)modal[row]+=M.val[row][q]*br[M.col[row][q]];modal[row]*=M.row_scale[row];}
 ref->derivatives.raw.along(2,ref->derivatives.inverse,4,modal.data(),physical.data());double inverse_residual=relative(balanced,physical);
 fprintf(stderr,"fields %.3e JVP %.3e M difference %.3e balanced M residual %.3e\n",fields_error,jvp_error,inverse_error,inverse_residual);
 CHECK(inverse_residual<=1e-9);
 HiSpID_Config small=c;small.memory_limit_mib=16;small.n[0]=80;small.n[1]=160;small.n[2]=16;CHECK(!HiSpID_create_with_execution(&small,1));
 // Find a budget accepted by the inherited CPU estimate but rejected by
 // the added backend estimate, before allocating any large context.
 small=c;bool backend_rejection=false;
 for(int mib=16;mib<=128;mib++){small.memory_limit_mib=mib;if(hispid::valid(small,false)){auto*rejected=HiSpID_create_with_execution(&small,1);if(!rejected){CHECK(std::string(HiSpID_last_error()).find("Kokkos aggregate")!=std::string::npos);backend_rejection=true;break;}HiSpID_destroy(rejected);}}
 CHECK(backend_rejection);
 printf("HiSpID %s actual kernels: %d checks; fields %.3e JVP %.3e M difference %.3e M residual %.3e\n",puncture::Exec::name(),checks,fields_error,jvp_error,inverse_error,inverse_residual);return 0;
}
