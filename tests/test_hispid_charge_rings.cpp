// Private-TU controls for cached meridional values and first derivatives.
#include "../src/HiSpID_solver.cpp"
#include <iostream>
int main(){
 HiSpID_Config cfg;HiSpID_default_config(&cfg);cfg.n[0]=8;cfg.n[1]=10;cfg.n[2]=12;
 const double axis[3]={.3,.4,std::sqrt(.75)},offset[3]={.7,-.2,.4};
 cfg.hole[0].mass=.6;cfg.hole[1].mass=.4;
 for(int i=0;i<3;i++){cfg.hole[0].center[i]=offset[i]+3*axis[i];cfg.hole[1].center[i]=offset[i]-3*axis[i];}
 const double spins[2][3]={{.05,-.04,.09},{.03,-.02,.05}},vel[2][3]={{.08,.02,-.03},{-.07,.01,.05}};
 for(int h=0;h<2;h++)for(int i=0;i<3;i++){cfg.hole[h].spin[i]=spins[h][i];cfg.hole[h].velocity[i]=vel[h][i];}
 HiSpID_Data*s=HiSpID_create_sampler(&cfg);if(!s)return 1;
 std::vector<double>values(s->ntotal);for(int p=0;p<s->ntotal;p++)values[p]=1e-5*std::sin(.23*p+.17)+2e-6*std::cos(.41*p);
 if(HiSpID_set_unknowns(s,values.data(),values.size()))return 1;
 double field_error=0,metric_error=0,gradient_error=0;
 for(double mu:{-.79,.13,.63}){
  const double radius=30,rho=radius*std::sqrt(1-mu*mu);auto ring=sphere_ring(*s,radius*mu,rho);
  for(int k=0;k<33;k++){
   const double phi=2*Pi*(k+.37)/33,x[3]={radius*mu,rho*std::cos(phi),rho*std::sin(phi)};
   Fields raw{},cached{};sample_fields(*s,x,raw);sample_fields(*s,x,cached,&ring);
   for(int v=0;v<4;v++)for(int d=0;d<4;d++)field_error=std::max(field_error,std::abs(raw[v][d]-cached[v][d])/(1+std::abs(raw[v][d])));
   double lab[3];for(int i=0;i<3;i++){lab[i]=s->origin[i];for(int j=0;j<3;j++)lab[i]+=s->frame[i][j]*x[j];}
   HiSpID_Point a,b;double da[27],db[27];
   if(sample_with_derivatives(s,1,lab,&a,da)||sample_with_derivatives(s,1,lab,&b,db,&ring,x))return 1;
   for(int i=0;i<9;i++){
    metric_error=std::max({metric_error,std::abs(a.gamma[i]-b.gamma[i])/(1+std::abs(a.gamma[i])),std::abs(a.Kij[i]-b.Kij[i])/(1+std::abs(a.Kij[i]))});
   }
   for(int i=0;i<27;i++)gradient_error=std::max(gradient_error,std::abs(da[i]-db[i])/(1+std::abs(da[i])));
  }
 }
 double before[7],after[7];if(HiSpID_charges(s,s->origin,30,24,48,before))return 1;
 // Change a monopole explicitly; every subsequent call must rebuild its ring.
 for(int p=0;p<cfg.n[0]*cfg.n[1];p++)values[4*p]+=1e-4;
 if(HiSpID_set_unknowns(s,values.data(),values.size())||HiSpID_charges(s,s->origin,30,24,48,after))return 1;
 const double change=std::abs(after[0]-before[0]);HiSpID_destroy(s);
 std::cout<<"ring scalar/vector sine/cosine/Nyquist value-gradient error "<<field_error
          <<", translated/rotated physical metric-K error "<<metric_error
          <<", lab metric-gradient error "<<gradient_error
          <<", changed-unknown energy response "<<change<<'\n';
 return field_error<1e-12&&metric_error<1e-11&&gradient_error<1e-11&&change>1e-8?0:1;
}
