#include "HiSpID_internal.hpp"
#include "HiSpID_axis.hpp"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <map>
#include <numeric>
#include <cstdio>
#include <gsl/gsl_linalg.h>
#ifndef HISPID_INFINITY_EQUILIBRATION
#define HISPID_INFINITY_EQUILIBRATION 0
#endif
static_assert(HISPID_INFINITY_EQUILIBRATION==0||HISPID_INFINITY_EQUILIBRATION==1,"infinity equilibration must be0or1");
extern "C" {
#include "TwoPunctures.h"
}
using hispid::Jet;
namespace {
const int hi[6]={1,1,1,2,2,3},hj[6]={1,2,3,2,3,3};
using Fields=std::array<std::array<double,10>,4>;
struct Cached {
 double psi,R,K,g,lapPsi,divM[3],gradK[3],M[9],inv[9];
 double lap[10],vec[3][3][10],L[9][3][4];
 double far_correction[10];
 double weight;
};
/* First-order tensor algebra for divergence of L. The scalar/vector fields
 * carry second derivatives, but background metric/connection need only
 * first derivatives here. This is also the explicit nonmetric-compatible
 * definition used for Eq27/28 inside attenuation zones. */
void L_and_div(const Jet metric[3][3],const Jet inv[3][3],const Jet C[3][3][3],
               const Fields&u,double L[3][3],double*out){
 double db[3][3]={},ddb[3][3][3]={},div=0,ddiv[3]={};
 double bu[3][3][3]={};
 for(int k=0;k<3;k++)for(int d=0;d<6;d++)bu[k][hi[d]-1][hj[d]-1]=bu[k][hj[d]-1][hi[d]-1]=u[k+1][4+d];
 for(int i=0;i<3;i++)for(int k=0;k<3;k++){
  db[i][k]=u[k+1][i+1];for(int l=0;l<3;l++)db[i][k]+=(double)C[k][i][l].v*u[l+1][0];
  for(int d=0;d<3;d++){
   ddb[d][i][k]=bu[k][i][d];for(int l=0;l<3;l++)ddb[d][i][k]+=(double)C[k][i][l].d[d+1]*u[l+1][0]+(double)C[k][i][l].v*u[l+1][d+1];
  }
  if(i==k){div+=db[i][k];for(int d=0;d<3;d++)ddiv[d]+=ddb[d][i][k];}
 }
 double dL[3][3][3]={};
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  L[i][j]=-2.0/3*(double)metric[i][j].v*div;
  for(int k=0;k<3;k++)L[i][j]+=(double)metric[i][k].v*db[j][k]+(double)metric[j][k].v*db[i][k];
  for(int d=0;d<3;d++){
   dL[d][i][j]=-2.0/3*((double)metric[i][j].d[d+1]*div+(double)metric[i][j].v*ddiv[d]);
   for(int k=0;k<3;k++)dL[d][i][j]+=(double)metric[i][k].d[d+1]*db[j][k]+(double)metric[i][k].v*ddb[d][j][k]
     +(double)metric[j][k].d[d+1]*db[i][k]+(double)metric[j][k].v*ddb[d][i][k];
  }
 }
 if(!out)return;
 double up[3][3]={},dup[3][3][3]={};
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)for(int l=0;l<3;l++){
  up[i][j]+=(double)(inv[i][k].v*inv[j][l].v)*L[k][l];
  for(int d=0;d<3;d++)dup[d][i][j]+=(double)(inv[i][k].v*inv[j][l].v)*dL[d][k][l]
    +(double)(inv[i][k].d[d+1]*inv[j][l].v+inv[i][k].v*inv[j][l].d[d+1])*L[k][l];
 }
 for(int i=0;i<3;i++){
  out[i]=0;for(int j=0;j<3;j++){
   out[i]+=dup[j][i][j];for(int k=0;k<3;k++)out[i]+=(double)C[i][j][k].v*up[k][j]+(double)C[j][j][k].v*up[i][k];
  }
 }
}
void cache(const hispid::Background&b,Cached&c){
 std::memset(&c,0,sizeof(c));c.psi=b.psi.v;c.R=b.R;c.K=b.K.v;c.g=b.g.v;c.lapPsi=b.lapPsi;
 c.far_correction[0]=b.far_correction.v;
 for(int d=1;d<4;d++)c.far_correction[d]=b.far_correction.d[d];
 for(int d=0;d<6;d++)c.far_correction[d+4]=b.far_correction.h[hi[d]][hj[d]];
 for(int i=0;i<3;i++){
  c.divM[i]=b.divM[i];for(int j=0;j<3;j++){
   c.gradK[i]+=(double)(b.inv[i][j].v*b.K.d[j+1]);
   c.M[3*i+j]=b.M[i][j].v;c.inv[3*i+j]=b.inv[i][j].v;
  }
 }
 for(int d=0;d<10;d++){
  Jet u;if(d==0)u.v=1;else if(d<4)u.d[d]=1;else u.h[hi[d-4]][hj[d-4]]=u.h[hj[d-4]][hi[d-4]]=1;
  c.lap[d]=hispid::laplacian(b.opinv,b.opC,u);
  for(int k=0;k<3;k++){
   Fields f{};f[k+1][d]=1;double L[3][3],div[3];
   L_and_div(b.opmetric,b.opinv,b.opC,f,L,div);
   for(int i=0;i<3;i++)c.vec[i][k][d]=div[i];
   if(d<4){L_and_div(b.metric,b.inv,b.C,f,L,nullptr);for(int i=0;i<3;i++)for(int j=0;j<3;j++)c.L[3*i+j][k][d]=L[i][j];}
  }
 }
}
double contraction(const Cached&c,const double*A,const double*B){
 double q=0;for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)for(int l=0;l<3;l++)
  q+=c.inv[3*i+k]*c.inv[3*j+l]*A[3*i+j]*B[3*k+l];
 return q;
}
void eval(const Cached&c,const Fields&u,double*out,const Fields*du=nullptr){
 double psi=c.psi+u[0][0],A[9];if(!(psi>0)){for(int k=0;k<4;k++)out[k]=NAN;return;}
 for(int a=0;a<9;a++){
  A[a]=c.M[a];for(int k=0;k<3;k++)for(int d=0;d<4;d++)A[a]+=c.L[a][k][d]*u[k+1][d];
 }
 const Fields&v=du?*du:u;for(int k=0;k<4;k++)out[k]=0;
 for(int d=0;d<10;d++)out[0]+=c.lap[d]*v[0][d];
 for(int i=0;i<3;i++)for(int k=0;k<3;k++)for(int d=0;d<10;d++)out[i+1]+=c.vec[i][k][d]*v[k+1][d];
 double A2=contraction(c,A,A);
 if(!du){
  out[0]+=c.g*(-psi*c.R/8-std::pow(psi,5)*c.K*c.K/12+A2/(8*std::pow(psi,7))+c.lapPsi);
  for(int k=0;k<3;k++)out[k+1]+=c.g*(c.divM[k]-2.0/3*std::pow(psi,6)*c.gradK[k]);
 }else{
  double DA[9]={};for(int a=0;a<9;a++)for(int k=0;k<3;k++)for(int d=0;d<4;d++)DA[a]+=c.L[a][k][d]*v[k+1][d];
  out[0]+=c.g*((-c.R/8-5*std::pow(psi,4)*c.K*c.K/12-7*A2/(8*std::pow(psi,8)))*v[0][0]+contraction(c,A,DA)/(4*std::pow(psi,7)));
  for(int k=0;k<3;k++)out[k+1]-=c.g*4*std::pow(psi,5)*v[0][0]*c.gradK[k];
 }
}
}
struct HiSpID_Data {
 HiSpID_Config config,local;
 double origin[3],frame[3][3],b;
 int npt,ntotal;
 std::vector<Cached> geometry;
 std::vector<double> values,coefficients;
 derivs *work=nullptr;
 hispid::AxisDerivatives derivatives;
 std::vector<Fields> basefields;
 HiSpID_Diagnostics diag{};
 bool coefficients_valid=false;
 bool sampler_only=false;
};
namespace {
int pindex(const HiSpID_Data&s,int i,int j,int k){return i+s.local.n[0]*(j+s.local.n[1]*k);}
double dot(const std::vector<double>&a,const std::vector<double>&b){double q=0;for(size_t i=0;i<a.size();i++)q+=a[i]*b[i];return q;}
double norm2v(const std::vector<double>&v){return std::sqrt(dot(v,v));}
double norminf(const std::vector<double>&v){double q=0;for(double a:v){if(!std::isfinite(a))return INFINITY;q=std::max(q,std::abs(a));}return q;}
void gather(derivs*w,int p,Fields&f){
 double *a[10]={w->d0,w->d1,w->d2,w->d3,w->d11,w->d12,w->d13,w->d22,w->d23,w->d33};
 for(int k=0;k<4;k++)for(int d=0;d<10;d++)f[k][d]=a[d][4*p+k];
}
void transform(HiSpID_Data&s,int i,int j,int k,Fields&f,double*x=nullptr){
 double A=s.derivatives.coordinate[0][i],B=s.derivatives.coordinate[1][j];
 double dat[10][4];derivs w{};double **ptr[10]={&w.d0,&w.d1,&w.d2,&w.d3,&w.d11,&w.d12,&w.d13,&w.d22,&w.d23,&w.d33};
 for(int d=0;d<10;d++){*ptr[d]=dat[d];for(int v=0;v<4;v++)dat[d][v]=(A-1)*f[v][d];}
 for(int v=0;v<4;v++){dat[1][v]+=f[v][0];dat[4][v]+=2*f[v][1];dat[5][v]+=f[v][2];dat[6][v]+=f[v][3];}
 double X,R,xx,rr,y,z;AB_To_XR(4,A,B,&X,&R,&w);C_To_c(4,X,R,&xx,&rr,s.b,&w);
 rx3_To_xyz(4,xx,rr,2*Pi*k/s.local.n[2],&y,&z,&w);
 for(int v=0;v<4;v++)for(int d=0;d<10;d++)f[v][d]=dat[d][v];
 if(x){x[0]=xx;x[1]=y;x[2]=z;}
}
void fields(HiSpID_Data&s,const double*values,std::vector<Fields>&f,bool include_reference=true){
 std::copy(values,values+s.ntotal,s.work->d0);
 s.derivatives.apply(4,s.work);
 f.resize(s.npt);for(int k=0;k<s.local.n[2];k++)for(int j=0;j<s.local.n[1];j++)for(int i=0;i<s.local.n[0];i++){
  int p=pindex(s,i,j,k);gather(s.work,p,f[p]);transform(s,i,j,k,f[p]);
  if(include_reference)for(int d=0;d<10;d++)f[p][0][d]+=s.geometry[p].far_correction[d];
 }
}
void residual(HiSpID_Data&s,const double*v,double*r){
 fields(s,v,s.basefields);for(int p=0;p<s.npt;p++){
  eval(s.geometry[p],s.basefields[p],r+4*p);for(int k=0;k<4;k++)r[4*p+k]*=s.geometry[p].weight;
 }
}
void jvp(HiSpID_Data&s,const double*v,double*r){
 std::vector<Fields>d;fields(s,v,d,false);for(int p=0;p<s.npt;p++){
  eval(s.geometry[p],s.basefields[p],r+4*p,&d[p]);for(int k=0;k<4;k++)r[4*p+k]*=s.geometry[p].weight;
 }
}
/* Exact block elimination of the five-point modal FD approximation. Each
 * polar row is a dense radial block after elimination. Cosine/sine partners
 * share factors, as do the three approximate vector equations. This removes
 * the long-wavelength error of ILU(0); residuals/JVPs remain pseudospectral. */
struct ModalBlock {
 int na=0,nb=0;
 std::vector<double>lu,transfer,lower,upper;
 std::vector<size_t>permutation;
 void factor(){
  std::vector<double>rhs(na);int sign=0;
  for(int j=0;j<nb;j++){
   double*block=lu.data()+j*na*na;
   if(j)for(int i=0;i<na;i++)for(int q=0;q<na;q++)
    block[i*na+q]-=lower[j*na+i]*transfer[((j-1)*na+i)*na+q];
   auto A=gsl_matrix_view_array(block,na,na);
   gsl_permutation p{size_t(na),permutation.data()+j*na};gsl_permutation_init(&p);
   if(gsl_linalg_LU_decomp(&A.matrix,&p,&sign))throw std::runtime_error("Modal block factorization failed");
   for(int i=0;i<na;i++)if(!std::isfinite(block[i*na+i])||std::abs(block[i*na+i])<1e-30)
    throw std::runtime_error("Modal block singular pivot");
   if(j+1<nb)for(int q=0;q<na;q++){
    std::fill(rhs.begin(),rhs.end(),0);rhs[q]=upper[j*na+q];
    auto b=gsl_vector_view_array(rhs.data(),na);
    auto x=gsl_vector_view_array_with_stride(transfer.data()+j*na*na+q,na,na);
    if(gsl_linalg_LU_solve(&A.matrix,&p,&b.vector,&x.vector))throw std::runtime_error("Modal transfer solve failed");
   }
  }
 }
 void solve(std::vector<double>&x,int mode,int component)const{
  std::vector<double>f(na*nb);
  for(int j=0;j<nb;j++){
   for(int i=0;i<na;i++)f[j*na+i]=x[4*(i+na*(j+nb*mode))+component]-(j?lower[j*na+i]*f[(j-1)*na+i]:0);
   auto A=gsl_matrix_const_view_array(lu.data()+j*na*na,na,na);
   gsl_permutation p{size_t(na),const_cast<size_t*>(permutation.data()+j*na)};
   auto y=gsl_vector_view_array(f.data()+j*na,na);
   if(gsl_linalg_LU_svx(&A.matrix,&p,&y.vector))throw std::runtime_error("Modal forward solve failed");
  }
  for(int j=nb-2;j>=0;j--)for(int i=0;i<na;i++){
   double sum=0;for(int q=0;q<na;q++)sum+=transfer[(j*na+i)*na+q]*f[(j+1)*na+q];f[j*na+i]-=sum;
  }
  for(int j=0;j<nb;j++)for(int i=0;i<na;i++)x[4*(i+na*(j+nb*mode))+component]=f[j*na+i];
 }
};
struct Sparse {
 std::vector<std::vector<int>>col;
 std::vector<std::vector<double>>val;
 std::vector<int>diag;
 const hispid::AxisDerivatives*modal=nullptr;
 std::vector<double>row_scale;
 std::vector<ModalBlock>blocks;
 int scalar_factorizations=0,vector_factorizations=0;
 void factor(std::vector<ModalBlock>*vector_cache=nullptr){
  scalar_factorizations=vector_factorizations=0;
  if(modal){
   const int na=modal->n[0],nb=modal->n[1],half=modal->n[2]/2;
   blocks.resize(2*(half+1));
   const bool reuse=vector_cache&&!vector_cache->empty();
   if(reuse&&vector_cache->size()!=size_t(half+1))throw std::runtime_error("Invalid modal vector cache");
   if(reuse)for(const auto&B:*vector_cache)
    if(B.na!=na||B.nb!=nb||B.lu.size()!=size_t(na*na*nb)||B.transfer.size()!=size_t(na*na*nb)||
       B.lower.size()!=size_t(na*nb)||B.upper.size()!=size_t(na*nb)||B.permutation.size()!=size_t(na*nb))
     throw std::runtime_error("Incomplete modal vector cache");
   for(int mode=0;mode<=half;mode++)for(int v=0;v<2;v++){
    auto&B=blocks[2*mode+v];B.na=na;B.nb=nb;
    if(v&&reuse){
     B=std::move((*vector_cache)[mode]);
     if(B.na!=na||B.nb!=nb)throw std::runtime_error("Modal vector cache grid mismatch");
     continue;
    }
    B.lu.assign(na*na*nb,0);B.transfer.assign(na*na*nb,0);
    B.lower.assign(na*nb,0);B.upper.assign(na*nb,0);B.permutation.resize(na*nb);
    for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
     const int row=4*(i+na*(j+nb*mode))+v;
     for(size_t q=0;q<col[row].size();q++){
      const int local=col[row][q]/4-na*nb*mode,ci=local%na,cj=local/na;
      if(cj==j)B.lu[(j*na+i)*na+ci]+=val[row][q];
      else if(ci==i&&cj==j-1)B.lower[j*na+i]+=val[row][q];
      else if(ci==i&&cj==j+1)B.upper[j*na+i]+=val[row][q];
      else throw std::runtime_error("Modal FD stencil is not block tridiagonal");
     }
    }
    B.factor();
    if(v)vector_factorizations++;else scalar_factorizations++;
   }
   if(vector_cache)vector_cache->clear();
   return;
  }
  int n=col.size();diag.resize(n);
  for(int i=0;i<n;i++){
   auto&cc=col[i];auto&vv=val[i];int di=std::lower_bound(cc.begin(),cc.end(),i)-cc.begin();diag[i]=di;
   for(int a=0;a<di;a++){
    int j=cc[a];double piv=val[j][diag[j]];if(std::abs(piv)<1e-30)throw std::runtime_error("ILU zero pivot");
    vv[a]/=piv;
    for(size_t q=diag[j]+1;q<col[j].size();q++){
     auto it=std::lower_bound(cc.begin()+a+1,cc.end(),col[j][q]);if(it!=cc.end()&&*it==col[j][q])vv[it-cc.begin()]-=vv[a]*val[j][q];
    }
   }
   if(!std::isfinite(vv[di])||std::abs(vv[di])<1e-30)throw std::runtime_error("ILU singular pivot");
  }
 }
 void retain_vectors(std::vector<ModalBlock>&cache){
  if(!modal)throw std::runtime_error("Vector reuse requires modal blocks");
  const int half=modal->n[2]/2;cache.resize(half+1);
  for(int mode=0;mode<=half;mode++)cache[mode]=std::move(blocks[2*mode+1]);
 }
 void solve(const std::vector<double>&b,std::vector<double>&x)const{
  int n=col.size();x=b;
  if(modal){
   const int N=modal->n[2],stride=modal->n[0]*modal->n[1];
   for(int line=0;line<stride;line++)for(int v=0;v<4;v++){
    double sum=0,error=0;for(int k=0;k<N;k++){double y=b[4*(line+k*stride)+v]-error,t=sum+y;error=(t-sum)-y;sum=t;}
    const double mean=sum/N;
    for(int mode=0;mode<N;mode++){
     double value=mode==0?sum/std::sqrt(double(N)):0;
     if(mode>0)for(int k=0;k<N;k++)value+=modal->forward[mode*N+k]*(b[4*(line+k*stride)+v]-mean);
     const int row=4*(line+mode*stride)+v;x[row]=value/row_scale[row];
    }
   }
   const int half=N/2;
   for(int k=0;k<N;k++)for(int v=0;v<4;v++)
    blocks[2*(k<=half?k:k-half)+(v?1:0)].solve(x,k,v);
   return;
  }
  for(int i=0;i<n;i++)for(int k=0;k<diag[i];k++)x[i]-=val[i][k]*x[col[i][k]];
  for(int i=n-1;i>=0;i--){for(size_t k=diag[i]+1;k<col[i].size();k++)x[i]-=val[i][k]*x[col[i][k]];x[i]/=val[i][diag[i]];}
 }
};
// Modal flat-Laplacian block preconditioner for the regular P unknowns.
// B_rm follows the analytic prolate Laplacian. Curved coefficients are
// represented by their azimuthal average; the vector block uses4/3 Delta.
struct AzimuthalAverage { double mu=0,potential=0; };
AzimuthalAverage azimuthal_average(const HiSpID_Data&s,int i,int j){
 AzimuthalAverage out;const int np=s.local.n[2];
 for(int phi=0;phi<np;phi++){
  int q=pindex(s,i,j,phi);const auto&g=s.geometry[q];const auto&u=s.basefields[q];double A[9];
  for(int x=0;x<9;x++){A[x]=g.M[x];for(int v=0;v<3;v++)for(int d=0;d<4;d++)A[x]+=g.L[x][v][d]*u[v+1][d];}
  const double psi=g.psi+u[0][0];out.mu+=(g.inv[0]+g.inv[4]+g.inv[8])/(3*np);
  out.potential+=g.g*(-g.R/8-5*std::pow(psi,4)*g.K*g.K/12-7*contraction(g,A,A)/(8*std::pow(psi,8)))/np;
 }
 return out;
}
Sparse preconditioner(HiSpID_Data&s,std::vector<ModalBlock>*vector_cache=nullptr,bool share_averages=true){
 Sparse mat;mat.col.resize(s.ntotal);mat.val.resize(s.ntotal);mat.row_scale.resize(s.ntotal);mat.modal=&s.derivatives;
 const int na=s.local.n[0],nb=s.local.n[1],np=s.local.n[2],half=np/2;
 const double ha=Pi/na,hb=Pi/nb;
 std::vector<AzimuthalAverage>averages;
 if(share_averages){
  averages.resize(na*nb);
  for(int j=0;j<nb;j++)for(int i=0;i<na;i++)averages[i+na*j]=azimuthal_average(s,i,j);
 }
 for(int k=0;k<np;k++)for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
  const int mode=k<=half?k:k-half,r=hispid::AxisDerivatives::exponent(mode),p=pindex(s,i,j,k);
  const double a=.5*(s.derivatives.coordinate[0][i]+1),t=a*a,B=s.derivatives.coordinate[1][j],eta=-2*B/(1+B*B),s2=1-eta*eta,h2=4*t/((1-t)*(1-t)),D=s.b*s.b*(h2+s2);
  const double c=-2*(1-t)*std::pow(a*std::sqrt(s2),r),weight=s.geometry[pindex(s,i,j,0)].weight;
  const auto average=share_averages?averages[i+na*j]:azimuthal_average(s,i,j);
  const double mu=average.mu,potential=average.potential;
  const double ctt=t*(1-t)*(1-t),ct=(1-t)*((r+1)*(1-t)-2*t),cee=s2,ce=-2*(r+1)*eta;
  const double c0=-(r+1)*(r+1-t)+(r*r-mode*mode)*(1/h2+1/s2);
  std::map<int,double>stencil;
  auto add=[&](int di,int dj,double v){int q=Index(0,i+di,j+dj,0,1,na,nb,1);stencil[q]+=v;};
  add(0,0,c0);
  const double al=ha*(i+.5),be=hb*(j+.5),sa=std::sin(al),sb=std::sin(be);
  const double lambda=hispid::AxisDerivatives::radial_stretch,sigma=.5*(1-std::cos(al)),d=1-(1-lambda)*sigma,dt=lambda/(d*d),ddt=2*(1-lambda)*lambda/(d*d*d),ta=.5*sa*dt,taa=.5*std::cos(al)*dt+.25*sa*sa*ddt;
  const double kappa=hispid::AxisDerivatives::angular_stretch,zeta=-std::cos(be),T=std::tanh(kappa),tanh=std::tanh(kappa*zeta),de=kappa*(1-tanh*tanh)/T,dde=-2*kappa*kappa*tanh*(1-tanh*tanh)/T,eb=de*sb,ebb=de*std::cos(be)+dde*sb*sb;
  const double aa=ctt/(ta*ta),ab=ct/ta-ctt*taa/(ta*ta*ta);
  const double ba=cee/(eb*eb),bb=ce/eb-cee*ebb/(eb*eb*eb);
  add(1,0,aa/(ha*ha)+ab/(2*ha));add(-1,0,aa/(ha*ha)-ab/(2*ha));add(0,0,-2*aa/(ha*ha));
  add(0,1,ba/(hb*hb)+bb/(2*hb));add(0,-1,ba/(hb*hb)-bb/(2*hb));add(0,0,-2*ba/(hb*hb));
  for(int v=0;v<4;v++){
   const int row=4*p+v;mat.row_scale[row]=weight*mu*(v?4./3:1.)*c/D;
   for(auto entry:stencil){int col=4*(entry.first+na*nb*k)+v;double value=entry.second;
    if(col==row&&v==0)value+=potential*D/mu;
    if(value!=0||col==row){mat.col[row].push_back(col);mat.val[row].push_back(value);}
   }
  }
 }
 mat.factor(vector_cache);return mat;
}
bool gmres(HiSpID_Data&s,const Sparse&M,const std::vector<double>&rhs,std::vector<double>&x,double rtol, bool eager=false){
 int n=s.ntotal,m=s.local.krylov_restart;x.assign(n,0);
 std::vector<double>Ax(n),r(n),w(n);double target=rtol*norm2v(rhs),beta=norm2v(rhs);
 if(beta<=target)return true;
 // Retain visited columns across restarts, including zero-initialized breakdown columns.
 // The eager path is a private equivalence control; the public solve always grows lazily.
 std::vector<std::vector<double>>V(m+1),Z(m);V[0].resize(n);
 if(eager){for(auto&v:V)v.resize(n);for(auto&z:Z)z.resize(n);}
 std::vector<std::vector<double>>H(m+1,std::vector<double>(m));std::vector<double>cs(m),sn(m),g(m+1),y(m);
 int total=0;while(total<s.local.max_krylov){
  jvp(s,x.data(),Ax.data());for(int a=0;a<n;a++)r[a]=rhs[a]-Ax[a];beta=norm2v(r);
  if(beta<=target)return true;if(!std::isfinite(beta))return false;
  for(int a=0;a<n;a++)V[0][a]=r[a]/beta;std::fill(g.begin(),g.end(),0);g[0]=beta;
  for(auto &h:H)std::fill(h.begin(),h.end(),0);
  int used=0;
  for(int k=0;k<m&&total<s.local.max_krylov;k++){
   if(Z[k].empty())Z[k].resize(n);if(V[k+1].empty())V[k+1].resize(n);
   M.solve(V[k],Z[k]);jvp(s,Z[k].data(),w.data());
   /* Two-pass modified Gram-Schmidt limits loss of orthogonality. */
   for(int pass=0;pass<2;pass++)for(int j=0;j<=k;j++){
    double a=dot(w,V[j]);H[j][k]+=a;for(int q=0;q<n;q++)w[q]-=a*V[j][q];
   }
   H[k+1][k]=norm2v(w);if(H[k+1][k]>0)for(int q=0;q<n;q++)V[k+1][q]=w[q]/H[k+1][k];
   for(int j=0;j<k;j++){double a=cs[j]*H[j][k]+sn[j]*H[j+1][k];H[j+1][k]=-sn[j]*H[j][k]+cs[j]*H[j+1][k];H[j][k]=a;}
   double dd=std::hypot(H[k][k],H[k+1][k]);if(!(dd>0))return false;
   cs[k]=H[k][k]/dd;sn[k]=H[k+1][k]/dd;H[k][k]=dd;H[k+1][k]=0;
   g[k+1]=-sn[k]*g[k];g[k]*=cs[k];used=k+1;total++;s.diag.krylov_iterations++;
   if(std::abs(g[k+1])<=target)break;
  }
  for(int j=used-1;j>=0;j--){y[j]=g[j];for(int k=j+1;k<used;k++)y[j]-=H[j][k]*y[k];y[j]/=H[j][j];}
  for(int j=0;j<used;j++)for(int q=0;q<n;q++)x[q]+=Z[j][q]*y[j];
 }
 jvp(s,x.data(),Ax.data());for(int a=0;a<n;a++)r[a]=rhs[a]-Ax[a];return norm2v(r)<=target;
}
void update_diag(HiSpID_Data&s,const std::vector<double>&res){
 for(int k=0;k<4;k++)s.diag.scaled_linf[k]=s.diag.unscaled_linf[k]=0;
 for(int p=0;p<s.npt;p++)for(int k=0;k<4;k++){
  if(!std::isfinite(res[4*p+k])){s.diag.scaled_linf[k]=s.diag.unscaled_linf[k]=INFINITY;continue;}
  s.diag.scaled_linf[k]=std::max(s.diag.scaled_linf[k],std::abs(res[4*p+k]));
  s.diag.unscaled_linf[k]=std::max(s.diag.unscaled_linf[k],std::abs(res[4*p+k]/s.geometry[p].weight));
 }
}
void to_local(const HiSpID_Data&s,const double*x,double*y){
 for(int i=0;i<3;i++){y[i]=0;for(int j=0;j<3;j++)y[i]+=s.frame[j][i]*(x[j]-s.origin[j]);}
}
void make_coefficients(HiSpID_Data&s){
 if(s.coefficients_valid)return;
 s.coefficients.resize(s.ntotal);std::vector<double>temporary(s.ntotal);
 s.derivatives.raw.along(0,s.derivatives.coefficient[0],4,s.values.data(),temporary.data());
 s.derivatives.raw.along(1,s.derivatives.coefficient[1],4,temporary.data(),s.coefficients.data());
 s.coefficients_valid=true;
}
void polynomial_basis(int N,double x,std::vector<double>&value,std::vector<double>&first){
 value.resize(N);first.resize(N);value[0]=1;first[0]=0;value[1]=x;first[1]=1;
 for(int i=2;i<N;i++){value[i]=2*x*value[i-1]-value[i-2];first[i]=2*value[i-1]+2*x*first[i-1]-first[i-2];}
 value[0]=.5;
}
// Return P, P_t, P_eta for one orthonormal Fourier mode.
void evaluate_mode(const HiSpID_Data&s,int k,double t,double eta,double P[4][3]){
 const int na=s.local.n[0],nb=s.local.n[1];std::vector<double>Ta,Da,Tb,Db;
 const double lambda=hispid::AxisDerivatives::radial_stretch,den=lambda+(1-lambda)*t,sigma=t/den,kappa=hispid::AxisDerivatives::angular_stretch,T=std::tanh(kappa),zeta=std::atanh(eta*T)/kappa;
 polynomial_basis(na,2*sigma-1,Ta,Da);polynomial_basis(nb,zeta,Tb,Db);std::memset(P,0,12*sizeof(double));
 for(int j=0;j<nb;j++)for(int i=0;i<na;i++)for(int v=0;v<4;v++){
  const double c=s.coefficients[v+4*(i+na*(j+nb*k))];P[v][0]+=c*Ta[i]*Tb[j];P[v][1]+=(2*lambda/(den*den))*c*Da[i]*Tb[j];P[v][2]+=(T/(kappa*(1-T*T*eta*eta)))*c*Ta[i]*Db[j];
 }
}
using ModeValues=std::array<std::array<double,3>,4>;
struct MeridionalValues {
 double sx2,sr2,sx,cx,sr,cr,a,B,rho;
 std::vector<ModeValues> modes;
};
// A sphere centered at the map origin and aligned to its axis has identical
// meridional coordinates around each polar ring. Only the Fourier sum and
// seed geometry vary with phi; the expensive two-dimensional P sums do not.
MeridionalValues sphere_ring(HiSpID_Data&s,double axial,double rho){
 make_coefficients(s);MeridionalValues out;out.rho=rho;
 const double w0=axial/s.b,q=rho/s.b,t=.5*(w0*w0+q*q-1),root=std::hypot(t,q);
 out.sx2=t<0?q*q/(root-t):root+t;out.sr2=t>0?q*q/(root+t):root-t;
 out.sx=std::sqrt(out.sx2);out.cx=std::sqrt(1+out.sx2);out.sr=std::sqrt(out.sr2);
 out.cr=w0/out.cx;out.a=out.sx/(out.cx+1);out.B=-out.cr/(1+out.sr);
 const double eta=-2*out.B/(1+out.B*out.B);out.modes.resize(s.local.n[2]);
 for(int k=0;k<s.local.n[2];k++){
  double values[4][3];evaluate_mode(s,k,out.a*out.a,eta,values);
  for(int v=0;v<4;v++)for(int d=0;d<3;d++)out.modes[k][v][d]=values[v][d];
 }
 return out;
}
void evaluate_coefficients(const HiSpID_Data&s,double a,double B,double phi,double V[4][4],double stable_sinR=-1,const MeridionalValues*ring=nullptr){
 const int np=s.local.n[2],half=np/2;const double t=a*a,eta=-2*B/(1+B*B),sn=stable_sinR>=0?stable_sinR:(1-B*B)/(1+B*B);
 std::memset(V,0,16*sizeof(double));
 for(int k=0;k<np;k++){
  const int m=k<=half?k:k-half,r=hispid::AxisDerivatives::exponent(m);double P[4][3];
  if(ring){for(int v=0;v<4;v++)for(int d=0;d<3;d++)P[v][d]=ring->modes[k][v][d];}
  else evaluate_mode(s,k,t,eta,P);
  const double norm=std::sqrt((m==0||m==half?1.:2.)/np),F=norm*(k<=half?std::cos(m*phi):std::sin(m*phi)),DF=norm*m*(k<=half?-std::sin(m*phi):std::cos(m*phi));
  const double S=(1+a)*std::pow(a*sn,r),Sa=.5*(std::pow(a,r)+ (r? r*(1+a)*std::pow(a,r-1):0))*std::pow(sn,r);
  const double etaB=-2*sn/(1+B*B),snB=-4*B/std::pow(1+B*B,2);
  const double Sb=r?(1+a)*std::pow(a,r)*r*std::pow(sn,r-1)*snB:0;
  for(int v=0;v<4;v++){
   V[v][0]+=S*P[v][0]*F;V[v][1]+=(Sa*P[v][0]+S*a*P[v][1])*F;
   V[v][2]+=(Sb*P[v][0]+S*etaB*P[v][2])*F;V[v][3]+=S*P[v][0]*DF;
  }
 }
}
void sample_fields(HiSpID_Data&s,const double*x,Fields&f,const MeridionalValues*ring=nullptr){
 make_coefficients(s);
 const double rho=ring?ring->rho:std::hypot(x[1],x[2]);
 /* Exact Cartesian first-derivative limits of the C2 modal interpolant.
  * Only m0 and m1 contribute. A tiny axis band avoids losing transverse
  * coordinates to rounding during coordinate inversion. */
 if(rho<1e-10*s.b){
  double A=-1,B,axial,transverse;const bool inner=std::abs(x[0])<s.b;
  if(inner){
   const double R=std::acos(std::clamp(x[0]/s.b,-1.0,1.0));
   B=std::tan(R/2-Piq);
   axial=-(1+B*B)/(2*s.b*std::sin(R));transverse=1/(s.b*std::sin(R));
  }else{
   if(std::abs(x[0])==s.b)throw std::runtime_error("puncture map focus is excluded from correction sampling");
   const double X=std::acosh(std::abs(x[0])/s.b);
   A=2*std::tanh(X/2)-1;B=x[0]>0?-1:1;
   const double a=.5*(A+1),sign=x[0]>0?1:-1;
   axial=sign*(1-a*a)/(s.b*std::sinh(X));transverse=sign/(s.b*std::sinh(X));
  }
  double cosine[4][4],sine[4][4];evaluate_coefficients(s,.5*(A+1),B,0,cosine);evaluate_coefficients(s,.5*(A+1),B,Pih,sine);
  for(int v=0;v<4;v++){
   f[v][0]=(A-1)*cosine[v][0];
   f[v][1]=inner?(A-1)*cosine[v][2]*axial:(cosine[v][0]+(A-1)*cosine[v][1])*axial;
   f[v][2]=inner?(cosine[v][0]+(A-1)*cosine[v][1])*transverse:(A-1)*cosine[v][2]*transverse;
   f[v][3]=inner?(sine[v][0]+(A-1)*sine[v][1])*transverse:(A-1)*sine[v][2]*transverse;
   f[v][0]+=f[v][2]*x[1]+f[v][3]*x[2];
  }
  return;
 }
 /* Stable inverse of x=b cosh(X) cos(R), rho=b sinh(X) sin(R).
  * Rationalize the small root rather than subtracting nearly equal
  * puncture distances (which loses R close to the axis). */
 const double w0=x[0]/s.b,q=rho/s.b,t=.5*(w0*w0+q*q-1),root=std::hypot(t,q);
 const double sx2=ring?ring->sx2:(t<0?q*q/(root-t):root+t);
 const double sr2=ring?ring->sr2:(t>0?q*q/(root+t):root-t);
 const double sx=ring?ring->sx:std::sqrt(sx2),cx=ring?ring->cx:std::sqrt(1+sx2),sr=ring?ring->sr:std::sqrt(sr2),cr=ring?ring->cr:w0/cx,a=ring?ring->a:sx/(cx+1),B=ring?ring->B:-cr/(1+sr);
 double phi=std::atan2(x[2],x[1]);if(phi<0)phi+=2*Pi;
 double V[4][4];evaluate_coefficients(s,a,B,phi,V,sr,ring);
 const double denominator=s.b*(sx2+sr2),AX=1-a*a,BR=.5*(1+B*B),co=x[1]/rho,si=x[2]/rho;
 for(int v=0;v<4;v++){
  const double uA=V[v][0]-2*(1-a)*V[v][1],uB=-2*(1-a)*V[v][2],up=-2*(1-a)*V[v][3];
  f[v][0]=-2*(1-a)*V[v][0];f[v][1]=(uA*AX*sx*cr-uB*BR*cx*sr)/denominator;
  const double radial=(uA*AX*cx*sr+uB*BR*sx*cr)/denominator;
  f[v][2]=radial*co-up*si/rho;f[v][3]=radial*si+up*co/rho;
 }

}
}
extern "C" {
const char *HiSpID_residual_scaling(){
 return HISPID_INFINITY_EQUILIBRATION?"sin6_alpha_beta_times_one_minus_t_pow_minus6":"sin6_alpha_beta";
}
const char *HiSpID_unknown_parameterization(){
 if constexpr(hispid::AxisDerivatives::radial_stretch==.2&&hispid::AxisDerivatives::angular_stretch==2.)return "modal_P_C2prolate_mapped_v2";
 static const auto identifier=[](){std::array<char,128> value{};
  std::snprintf(value.data(),value.size(),"modal_P_C2prolate_map_v3_r%.17g_k%.17g",hispid::AxisDerivatives::radial_stretch,hispid::AxisDerivatives::angular_stretch);return value;}();
 return identifier.data();
}
int HiSpID_collocation_maps(double*out){
 if(!out)return -1;out[0]=hispid::AxisDerivatives::radial_stretch;out[1]=hispid::AxisDerivatives::angular_stretch;return 0;
}
static HiSpID_Data *create_context(const HiSpID_Config*c,bool sampler_only){
 hispid::last_error.clear();
 if(!c||!hispid::valid(*c,sampler_only)){hispid::last_error="invalid HiSpID configuration";return nullptr;}
 HiSpID_Data*s=nullptr;try{
  s=new HiSpID_Data;s->config=s->local=*c;s->sampler_only=sampler_only;
  double sep[3],len=0;for(int i=0;i<3;i++){s->origin[i]=.5*(c->hole[0].center[i]+c->hole[1].center[i]);sep[i]=c->hole[0].center[i]-c->hole[1].center[i];len+=sep[i]*sep[i];}
  len=std::sqrt(len);
  if(!std::isfinite(len))throw std::runtime_error("nonfinite puncture separation");
  if(len<1e-12)throw std::runtime_error("puncture map centers must be separated, including an inactive focus");
  s->b=len/2;for(int i=0;i<3;i++)s->frame[i][0]=sep[i]/len;
  int least=0;for(int i=1;i<3;i++)if(std::abs(s->frame[i][0])<std::abs(s->frame[least][0]))least=i;
  double norm=std::sqrt(1-s->frame[least][0]*s->frame[least][0]);
  for(int i=0;i<3;i++)s->frame[i][1]=((i==least?1.0:0.0)-s->frame[i][0]*s->frame[least][0])/norm;
  for(int i=0;i<3;i++)s->frame[i][2]=s->frame[(i+1)%3][0]*s->frame[(i+2)%3][1]-s->frame[(i+2)%3][0]*s->frame[(i+1)%3][1];
  for(int h=0;h<2;h++){
   to_local(*s,c->hole[h].center,s->local.hole[h].center);
   for(int i=0;i<3;i++){s->local.hole[h].spin[i]=s->local.hole[h].velocity[i]=0;for(int j=0;j<3;j++){
    s->local.hole[h].spin[i]+=s->frame[j][i]*c->hole[h].spin[j];s->local.hole[h].velocity[i]+=s->frame[j][i]*c->hole[h].velocity[j];}}
  }
  s->npt=c->n[0]*c->n[1]*c->n[2];s->ntotal=4*s->npt;s->values.resize(s->ntotal);
  s->diag.npoints=s->npt;
  if(sampler_only){s->derivatives.initialize(c->n,true);return s;}
  s->geometry.resize(s->npt);allocate_derivs(&s->work,s->ntotal);
  s->derivatives.initialize(c->n);
  for(int k=0;k<c->n[2];k++)for(int j=0;j<c->n[1];j++)for(int i=0;i<c->n[0];i++){
   Fields f{};double x[3];transform(*s,i,j,k,f,x);hispid::Background bg;hispid::background(s->local,x,bg);
   int p=pindex(*s,i,j,k);cache(bg,s->geometry[p]);
   double sn=std::sin(Pih*(2*i+1)/c->n[0])*std::sin(Pih*(2*j+1)/c->n[1]);s->geometry[p].weight=std::pow(sn,6);
   if constexpr(HISPID_INFINITY_EQUILIBRATION){
    const double a=.5*(s->derivatives.coordinate[0][i]+1);
    s->geometry[p].weight/=std::pow(1-a*a,6);
   }
  }
  s->diag.npoints=s->npt;
 }catch(const std::exception&e){hispid::last_error=e.what();HiSpID_destroy(s);return nullptr;}return s;
}
HiSpID_Data *HiSpID_create(const HiSpID_Config*c){return create_context(c,false);}
HiSpID_Data *HiSpID_create_sampler(const HiSpID_Config*c){return create_context(c,true);}
static bool sampling_context(HiSpID_Data*s){
 if(s&&s->sampler_only){hispid::last_error="sampling-only context cannot evaluate or solve collocation equations";return true;}
 return false;
}
int HiSpID_solve(HiSpID_Data*s){
 if(sampling_context(s))return -1;
 if(!s)return -1;hispid::last_error.clear();auto start=std::chrono::steady_clock::now();
 s->diag={};s->diag.npoints=s->npt;s->coefficients_valid=false;
 std::vector<double>r(s->ntotal),rhs(s->ntotal),step,trial(s->ntotal),rt(s->ntotal);
 // The normalized vector FD matrices depend only on this solve's fixed
 // map/grid/mode. Move their factors between Newton steps; scalar factors
 // are rebuilt from the current nonlinear potential at every step.
 std::vector<ModalBlock>vector_cache;
 try{
  residual(*s,s->values.data(),r.data());
  for(int it=0;it<=s->local.max_newton;it++){
   update_diag(*s,r);double err=norminf(r);if(err<=s->local.tolerance){s->diag.converged=1;break;}
   if(it==s->local.max_newton)break;
   s->diag.newton_iterations++;
   Sparse M=preconditioner(*s,&vector_cache);for(int i=0;i<s->ntotal;i++)rhs[i]=-r[i];
   double forcing=std::min(.05,std::max(1e-5,std::sqrt(err)));
   const bool linear_ok=gmres(*s,M,rhs,step,forcing);M.retain_vectors(vector_cache);
   if(!linear_ok){hispid::last_error="Krylov iteration limit";break;}
   bool accepted=false;double old=norm2v(r);
   for(double damping=1;damping>=1.0/1024;damping*=.5){
    for(int i=0;i<s->ntotal;i++)trial[i]=s->values[i]+damping*step[i];
    residual(*s,trial.data(),rt.data());
    if(norm2v(rt)<old*(1-1e-4*damping)&&std::isfinite(norminf(rt))){s->values=trial;r=rt;accepted=true;break;}
   }
   if(!accepted){hispid::last_error="Newton line search failed";residual(*s,s->values.data(),r.data());break;}
  }
 }catch(const std::exception&e){hispid::last_error=e.what();}
 residual(*s,s->values.data(),r.data());update_diag(*s,r);
 s->diag.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 if(!s->diag.converged&&hispid::last_error.empty())hispid::last_error="Newton iteration limit";
 return s->diag.converged?0:1;
}
int HiSpID_diagnostics(const HiSpID_Data*s,HiSpID_Diagnostics*d){if(!s||!d)return -1;*d=s->diag;return 0;}
int HiSpID_get_unknowns(const HiSpID_Data*s,double*v,int n){if(!s||!v||n!=s->ntotal)return -1;std::copy(s->values.begin(),s->values.end(),v);return 0;}
int HiSpID_set_unknowns(HiSpID_Data*s,const double*v,int n){if(!s||!v||n!=s->ntotal)return -1;for(int i=0;i<n;i++)if(!std::isfinite(v[i]))return -1;std::copy(v,v+n,s->values.begin());s->coefficients_valid=false;s->diag.converged=0;return 0;}
int HiSpID_residual(HiSpID_Data*s,const double*v,double*r){
 if(sampling_context(s))return -1;
 if(!s||!v||!r)return -1;
 try{residual(*s,v,r);}catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
int HiSpID_jvp(HiSpID_Data*s,const double*v,const double*d,double*r){
 if(sampling_context(s))return -1;
 if(!s||!v||!d||!r)return -1;
 try{fields(*s,v,s->basefields);jvp(*s,d,r);}catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
int HiSpID_equation_samples(HiSpID_Data*s,double*xyz,double*g,double*psi,double*HM){
 if(sampling_context(s))return -1;
 if(!s||!xyz||!g||!psi||!HM)return -1;
 try{
  fields(*s,s->values.data(),s->basefields);
  for(int k=0;k<s->local.n[2];k++)for(int j=0;j<s->local.n[1];j++)for(int i=0;i<s->local.n[0];i++){
   const int p=pindex(*s,i,j,k);const Cached&c=s->geometry[p];
   psi[p]=c.psi+s->basefields[p][0][0];g[p]=c.g;
   eval(c,s->basefields[p],HM+4*p);
   HM[4*p]*=-8/std::pow(psi[p],5);
   for(int d=1;d<4;d++)HM[4*p+d]/=std::pow(psi[p],10);
   double momentum[3]={};for(int a=0;a<3;a++)for(int b=0;b<3;b++)momentum[a]+=s->frame[a][b]*HM[4*p+b+1];
   for(int a=0;a<3;a++)HM[4*p+a+1]=momentum[a];
   Fields zero{};double x[3];transform(*s,i,j,k,zero,x);
   for(int a=0;a<3;a++){
    xyz[3*p+a]=s->origin[a];for(int b=0;b<3;b++)xyz[3*p+a]+=s->frame[a][b]*x[b];
   }
  }
 }catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
static int sample_with_derivatives(HiSpID_Data*s,int count,const double*xyz,HiSpID_Point*out,double*dgamma,const MeridionalValues*ring=nullptr,const double*local_xyz=nullptr){
 if(!s||!xyz||!out||count<0)return -1;
 try{for(int p=0;p<count;p++){
  double x[3];if(local_xyz)std::copy(local_xyz+3*p,local_xyz+3*p+3,x);else to_local(*s,xyz+3*p,x);
  hispid::Background bg;hispid::background(s->local,x,bg);Fields f{};sample_fields(*s,x,f,ring);
  f[0][0]+=bg.far_correction.v;
  double psi=bg.psi.v+f[0][0];if(!(psi>0))throw std::runtime_error("nonpositive solved conformal factor");
  double L[3][3],gam[3][3],K[3][3],A[3][3];L_and_div(bg.metric,bg.inv,bg.C,f,L,nullptr);
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   A[i][j]=bg.M[i][j].v+L[i][j];gam[i][j]=std::pow(psi,4)*(double)bg.metric[i][j].v;
   K[i][j]=A[i][j]/(psi*psi)+gam[i][j]*(double)bg.K.v/3;
  }
  std::memset(out+p,0,sizeof(*out));out[p].psi=psi;out[p].mean_curvature=bg.K.v;out[p].attenuation=bg.g.v;
  out[p].correction[0]=f[0][0];for(int i=0;i<3;i++)for(int a=0;a<3;a++)out[p].correction[i+1]+=s->frame[i][a]*f[a+1][0];
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int a=0;a<3;a++)for(int b=0;b<3;b++){
   double q=s->frame[i][a]*s->frame[j][b];out[p].gamma[3*i+j]+=q*gam[a][b];out[p].Kij[3*i+j]+=q*K[a][b];
   out[p].conformal_metric[3*i+j]+=q*(double)bg.metric[a][b].v;out[p].Atilde[3*i+j]+=q*A[a][b];
  }
  if(dgamma){
   std::fill(dgamma+27*p,dgamma+27*(p+1),0.0);
   for(int d=0;d<3;d++)for(int i=0;i<3;i++)for(int j=0;j<3;j++)
    for(int c=0;c<3;c++)for(int a=0;a<3;a++)for(int b=0;b<3;b++){
     const double dpsi=bg.psi.d[c+1]+f[0][c+1]+bg.far_correction.d[c+1];
     const double grad=4*std::pow(psi,3)*dpsi*(double)bg.metric[a][b].v
                        +std::pow(psi,4)*(double)bg.metric[a][b].d[c+1];
     dgamma[27*p+9*d+3*i+j]+=s->frame[d][c]*s->frame[i][a]*s->frame[j][b]*grad;
    }
  }
 }}catch(const std::exception&e){hispid::last_error=e.what();return -2;}return 0;
}
int HiSpID_sample(HiSpID_Data*s,int count,const double*xyz,HiSpID_Point*out){
 return sample_with_derivatives(s,count,xyz,out,nullptr);
}
int HiSpID_sample_with_derivatives(HiSpID_Data*s,int count,const double*xyz,HiSpID_Point*out,double*dgamma){
 if(!dgamma)return -1;
 return sample_with_derivatives(s,count,xyz,out,dgamma);
}
int HiSpID_charges(HiSpID_Data*s,const double*center,double radius,int nt,int np,double*out){
 if(!s||!center||!out||radius<=0||nt<4||np<8)return -1;
 try{
 std::fill(out,out+7,0.0);
 const bool centered=center[0]==s->origin[0]&&center[1]==s->origin[1]&&center[2]==s->origin[2];
 /* Gauss-Legendre cos(theta) and uniform phi. */
 for(int a=0;a<nt;a++){
  double mu=std::cos(Pi*(a+.75)/(nt+.5)),pp=0;
  for(int it=0;it<20;it++){
   double p0=1,p1=mu;for(int l=2;l<=nt;l++){double p=((2*l-1)*mu*p1-(l-1)*p0)/l;p0=p1;p1=p;}
   pp=nt*(mu*p1-p0)/(mu*mu-1);double step=p1/pp;mu-=step;if(std::abs(step)<1e-15)break;
  }double weight=2/((1-mu*mu)*pp*pp)*(2*Pi/np)*radius*radius;
  MeridionalValues ring;if(centered)ring=sphere_ring(*s,radius*mu,radius*std::sqrt(1-mu*mu));
  for(int b=0;b<np;b++){
   // Align the integration polar axis with the prolate separation axis.
   // This keeps high meridional polynomial degrees out of the azimuthal
   // quadrature. Normals/tensors and returned E,P,J remain in the lab frame.
   double ph=2*Pi*(b+.5)/np,local_normal[3]={mu,std::sqrt(1-mu*mu)*std::cos(ph),std::sqrt(1-mu*mu)*std::sin(ph)},n[3]={};
   for(int i=0;i<3;i++)for(int j=0;j<3;j++)n[i]+=s->frame[i][j]*local_normal[j];
   double x[3];for(int i=0;i<3;i++)x[i]=center[i]+radius*n[i];HiSpID_Point p;double dg[3][9];
   // Reuse the tested analytic physical metric gradient. This evaluates
   // the retained modal polynomial once rather than at13 FD stencil points.
   // Independent physical-FD charge comparisons remain validation controls.
   double local_x[3];for(int i=0;i<3;i++)local_x[i]=radius*local_normal[i];
   if(sample_with_derivatives(s,1,x,&p,&dg[0][0],centered?&ring:nullptr,centered?local_x:nullptr))return -2;
   double E=0,P[3]={};for(int i=0;i<3;i++)for(int j=0;j<3;j++){
    E+=n[i]*(dg[j][3*i+j]-dg[i][3*j+j]);
    P[i]+=(p.Kij[3*i+j]-p.mean_curvature*p.gamma[3*i+j])*n[j];
   }
   out[0]+=weight*E/(16*Pi);for(int i=0;i<3;i++){
    out[i+1]+=weight*P[i]/(8*Pi);out[i+4]+=weight*radius*(n[(i+1)%3]*P[(i+2)%3]-n[(i+2)%3]*P[(i+1)%3])/(8*Pi);
   }
  }
 }return 0;
 }catch(const std::exception&e){hispid::last_error=e.what();return -2;}
}
void HiSpID_destroy(HiSpID_Data*s){if(!s)return;if(s->work)free_derivs(s->work);delete s;}
}
