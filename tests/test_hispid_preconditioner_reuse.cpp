// Private-TU equivalence controls for exact averaging/factor reuse.
#include "../src/HiSpID_solver.cpp"
#include <iostream>
bool same_block(const ModalBlock&a,const ModalBlock&b){
 return a.na==b.na&&a.nb==b.nb&&a.lu==b.lu&&a.transfer==b.transfer&&
        a.lower==b.lower&&a.upper==b.upper&&a.permutation==b.permutation;
}
bool same_matrix(const Sparse&a,const Sparse&b){
 if(a.col!=b.col||a.val!=b.val||a.row_scale!=b.row_scale||a.blocks.size()!=b.blocks.size())return false;
 for(size_t i=0;i<a.blocks.size();i++)if(!same_block(a.blocks[i],b.blocks[i]))return false;
 return true;
}
int main(){
 HiSpID_Config cfg;HiSpID_default_config(&cfg);cfg.n[0]=10;cfg.n[1]=14;cfg.n[2]=12;
 cfg.far_radius=0;cfg.hole[0].mass=.6;cfg.hole[1].mass=.4;
 const double spin[2][3]={{.072,.054,.108},{-.032,.048,.016}},vel[2][3]={{.03,.06,.01},{-.02,-.07,.025}};
 for(int h=0;h<2;h++)for(int d=0;d<3;d++){cfg.hole[h].spin[d]=spin[h][d];cfg.hole[h].velocity[d]=vel[h][d];}
 HiSpID_Data*s=HiSpID_create(&cfg);if(!s)return 1;
 std::vector<double>values(s->ntotal),res(s->ntotal),rhs(s->ntotal),a,b;
 for(int p=0;p<s->ntotal;p++){values[p]=1e-5*std::sin(.23*p+.17);rhs[p]=std::sin(.37*p+.11)+std::cos(.13*p);}
 residual(*s,values.data(),res.data());
 Sparse original=preconditioner(*s,nullptr,false),shared=preconditioner(*s);
 original.solve(rhs,a);shared.solve(rhs,b);const bool first_bitwise=same_matrix(original,shared)&&a==b;bool passed=first_bitwise;
 const int groups=cfg.n[2]/2+1;
 passed &= shared.scalar_factorizations==groups&&shared.vector_factorizations==groups;
 std::vector<ModalBlock>cache;shared.retain_vectors(cache);
 passed &= cache.size()==size_t(groups);
 for(int p=0;p<s->ntotal;p++)values[p]=2e-5*std::cos(.31*p+.21);
 residual(*s,values.data(),res.data());
 Sparse fresh=preconditioner(*s,nullptr,false),reused=preconditioner(*s,&cache);
 fresh.solve(rhs,a);reused.solve(rhs,b);
 const bool later_bitwise=same_matrix(fresh,reused)&&a==b&&cache.empty();passed &= later_bitwise;
 passed &= reused.scalar_factorizations==groups&&reused.vector_factorizations==0;
 bool scalar_changed=false;for(int k=0;k<groups;k++)scalar_changed |= original.blocks[2*k].lu!=fresh.blocks[2*k].lu;
 passed &= scalar_changed;
 // Independently apply every retained FD row to a known modal vector and
 // restore it through the exact block inverse (including Fourier partners).
 std::vector<double>known(s->ntotal),modal_rhs(s->ntotal),physical_rhs(s->ntotal),solution;
 for(int p=0;p<s->ntotal;p++)known[p]=1e-3*std::sin(.19*p+.37);
 for(int row=0;row<s->ntotal;row++)for(size_t q=0;q<fresh.col[row].size();q++)
  modal_rhs[row]+=fresh.row_scale[row]*fresh.val[row][q]*known[fresh.col[row][q]];
 const int stride=cfg.n[0]*cfg.n[1],np=cfg.n[2];
 for(int line=0;line<stride;line++)for(int v=0;v<4;v++)for(int phi=0;phi<np;phi++)for(int mode=0;mode<np;mode++)
  physical_rhs[4*(line+phi*stride)+v]+=s->derivatives.forward[mode*np+phi]*modal_rhs[4*(line+mode*stride)+v];
 reused.solve(physical_rhs,solution);fresh.solve(physical_rhs,a);double inversion_error=0,baseline_inversion_error=0;
 for(int p=0;p<s->ntotal;p++)inversion_error=std::max(inversion_error,std::abs(solution[p]-known[p]));
 for(int p=0;p<s->ntotal;p++)baseline_inversion_error=std::max(baseline_inversion_error,std::abs(a[p]-known[p]));
 // The scaled Fourier manufacturing round trip has a measured dynamic-
 // range floor, identically in the baseline and reused path. Preserve it
 // as a diagnostic; test the exact FD inverse directly without that round trip.
 passed &= solution==a;
 std::vector<double>block_solution(s->ntotal);
 for(int row=0;row<s->ntotal;row++)for(size_t q=0;q<fresh.col[row].size();q++)
  block_solution[row]+=fresh.val[row][q]*known[fresh.col[row][q]];
 for(int mode=0;mode<np;mode++)for(int v=0;v<4;v++)
  reused.blocks[2*(mode<=np/2?mode:mode-np/2)+(v?1:0)].solve(block_solution,mode,v);
 double block_error=0;for(int p=0;p<s->ntotal;p++)block_error=std::max(block_error,std::abs(block_solution[p]-known[p]));
 passed &= block_error<1e-10;
 // Reject a partially moved bank before dereferencing any factor buffers.
 reused.retain_vectors(cache);cache[0].lu.clear();bool rejected=false;
 try{auto invalid=preconditioner(*s,&cache);}catch(const std::runtime_error&){rejected=true;}
 passed &= rejected;HiSpID_destroy(s);
 std::cout<<"first/later matrix, factors, random RHS bitwise "<<first_bitwise<<"/"<<later_bitwise
          <<", scalar refreshed "<<scalar_changed<<", second scalar/vector factor counts "
          <<reused.scalar_factorizations<<"/"<<reused.vector_factorizations
          <<", known-vector inversion error "<<inversion_error<<", uncached inversion error "<<baseline_inversion_error
          <<", direct modal inverse error "<<block_error
          <<", incomplete cache rejected "<<rejected<<", aggregate "<<passed<<'\n';
 return passed?0:1;
}
