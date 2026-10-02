#include "HiSpID_internal.hpp"
#include "HiSpID_spectral.hpp"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <map>
#include <numeric>
#include <cstdio>
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
 hispid::SpectralDerivatives derivatives;
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
 double A=-std::cos(Pih*(2*i+1)/s.local.n[0]),B=-std::cos(Pih*(2*j+1)/s.local.n[1]);
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
/* A row-local second-order FD approximation in uniform (alpha,beta,phi),
 * with reflected Chebyshev boundaries and periodic phi, used only as an
 * ILU(0) preconditioner. The actual solve/Jacobian are pseudospectral. */
struct Sparse {
 std::vector<std::vector<int>>col;
 std::vector<std::vector<double>>val;
 std::vector<int>diag;
 void factor(){
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
 void solve(const std::vector<double>&b,std::vector<double>&x)const{
  int n=col.size();x=b;
  for(int i=0;i<n;i++)for(int k=0;k<diag[i];k++)x[i]-=val[i][k]*x[col[i][k]];
  for(int i=n-1;i>=0;i--){for(size_t k=diag[i]+1;k<col[i].size();k++)x[i]-=val[i][k]*x[col[i][k]];x[i]/=val[i][diag[i]];}
 }
};
Sparse preconditioner(HiSpID_Data&s){
 Sparse mat;mat.col.resize(s.ntotal);mat.val.resize(s.ntotal);
 double step[3]={Pi/s.local.n[0],Pi/s.local.n[1],2*Pi/s.local.n[2]};
 for(int k=0;k<s.local.n[2];k++)for(int j=0;j<s.local.n[1];j++)for(int i=0;i<s.local.n[0];i++){
  int p=pindex(s,i,j,k);double al=step[0]*(i+.5),be=step[1]*(j+.5);
  double sn[2]={std::sin(al),std::sin(be)},co[2]={std::cos(al),std::cos(be)};
  std::map<int,std::array<double,10>>stencil;
  auto add=[&](int di,int dj,int dk,int d,double v){
   int q=Index(0,i+di,j+dj,k+dk,1,s.local.n[0],s.local.n[1],s.local.n[2]);stencil[q][d]+=v;
  };
  add(0,0,0,0,1);
  for(int dim=0;dim<3;dim++){
   int pos[3]={},neg[3]={};pos[dim]=1;neg[dim]=-1;
   add(pos[0],pos[1],pos[2],dim+1,.5/step[dim]);add(neg[0],neg[1],neg[2],dim+1,-.5/step[dim]);
   int d=dim==0?4:dim==1?7:9;
   add(pos[0],pos[1],pos[2],d,1/(step[dim]*step[dim]));add(neg[0],neg[1],neg[2],d,1/(step[dim]*step[dim]));add(0,0,0,d,-2/(step[dim]*step[dim]));
  }
  int mixed[3][3]={{0,1,5},{0,2,6},{1,2,8}};
  for(auto &mix:mixed)for(int sg1:{-1,1})for(int sg2:{-1,1}){
   int d[3]={};d[mix[0]]=sg1;d[mix[1]]=sg2;add(d[0],d[1],d[2],mix[2],sg1*sg2*.25/(step[mix[0]]*step[mix[1]]));
  }
  std::map<int,std::array<double,4>>rows[4];
  for(auto &entry:stencil){
   auto a=entry.second;
   a[4]=a[4]/(sn[0]*sn[0])-co[0]*a[1]/std::pow(sn[0],3);
   a[7]=a[7]/(sn[1]*sn[1])-co[1]*a[2]/std::pow(sn[1],3);
   a[5]/=sn[0]*sn[1];a[6]/=sn[0];a[8]/=sn[1];a[1]/=sn[0];a[2]/=sn[1];
   for(int v=0;v<4;v++){
    Fields d{};d[v]=a;transform(s,i,j,k,d);double out[4];eval(s.geometry[p],s.basefields[p],out,&d);
    for(int r=0;r<4;r++)rows[r][4*entry.first+v][0]+=out[r]*s.geometry[p].weight;
   }
  }
  for(int r=0;r<4;r++){
   int row=4*p+r;rows[r][row][0]+=0;
   for(auto &a:rows[r])if(a.second[0]!=0||a.first==row){mat.col[row].push_back(a.first);mat.val[row].push_back(a.second[0]);}
  }
 }
 mat.factor();return mat;
}
bool gmres(HiSpID_Data&s,const Sparse&M,const std::vector<double>&rhs,std::vector<double>&x,double rtol){
 int n=s.ntotal,m=s.local.krylov_restart;x.assign(n,0);
 std::vector<double>Ax(n),r(n),z(n),w(n);double target=rtol*norm2v(rhs),beta=norm2v(rhs);
 if(beta<=target)return true;
 std::vector<std::vector<double>>V(m+1,std::vector<double>(n)),Z(m,std::vector<double>(n));
 std::vector<std::vector<double>>H(m+1,std::vector<double>(m));std::vector<double>cs(m),sn(m),g(m+1),y(m);
 int total=0;while(total<s.local.max_krylov){
  jvp(s,x.data(),Ax.data());for(int a=0;a<n;a++)r[a]=rhs[a]-Ax[a];beta=norm2v(r);
  if(beta<=target)return true;if(!std::isfinite(beta))return false;
  for(int a=0;a<n;a++)V[0][a]=r[a]/beta;std::fill(g.begin(),g.end(),0);g[0]=beta;
  for(auto &h:H)std::fill(h.begin(),h.end(),0);
  int used=0;
  for(int k=0;k<m&&total<s.local.max_krylov;k++){
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
 s.coefficients.resize(s.ntotal);
 SpecCoef(s.local.n[0],s.local.n[1],s.local.n[2],4,s.values.data(),s.coefficients.data());
 s.coefficients_valid=true;
}
void chebyshev_basis(int n,double x,std::vector<double>&T,std::vector<double>&D){
 T.resize(n);D.resize(n);T[0]=1;D[0]=0;T[1]=x;D[1]=1;
 for(int i=2;i<n;i++){
  T[i]=2*x*T[i-1]-T[i-2];D[i]=2*T[i-1]+2*x*D[i-1]-D[i-2];
 }
 T[0]=.5;
}
void evaluate_coefficients(const HiSpID_Data&s,double A,double B,double phi,double V[4][4]){
 const int na=s.local.n[0],nb=s.local.n[1],np=s.local.n[2],m=np/2;
 std::vector<double>Ta,Da,Tb,Db;chebyshev_basis(na,A,Ta,Da);chebyshev_basis(nb,B,Tb,Db);
 std::memset(V,0,16*sizeof(double));
 for(int k=0;k<np;k++){
  const int frequency=k<=m?k:k-m;
  double F=k<=m?std::cos(frequency*phi):std::sin(frequency*phi);
  double DF=k<=m?-frequency*std::sin(frequency*phi):frequency*std::cos(frequency*phi);
  if(k==0||k==m){F*=.5;DF*=.5;}
  double vB[4][3]={};
  for(int j=0;j<nb;j++){
   double vA[4][2]={};
   for(int i=0;i<na;i++)for(int v=0;v<4;v++){
    const double c=s.coefficients[v+4*(i+na*(j+nb*k))];
    vA[v][0]+=c*Ta[i];vA[v][1]+=c*Da[i];
   }
   for(int v=0;v<4;v++){
    vB[v][0]+=vA[v][0]*Tb[j];vB[v][1]+=vA[v][1]*Tb[j];vB[v][2]+=vA[v][0]*Db[j];
   }
  }
  for(int v=0;v<4;v++){
   for(int d=0;d<3;d++)V[v][d]+=vB[v][d]*F;
   V[v][3]+=vB[v][0]*DF;
  }
 }
}
void sample_fields(HiSpID_Data&s,const double*x,Fields&f){
 make_coefficients(s);
 const double rho=std::hypot(x[1],x[2]);
 /* The Cartesian limit of the map-axis derivative is evaluated as in
  * thesis Sec.2.4.2, using four transverse points. Test the physical rho
  * before coordinate inversion so distant axis points cannot recurse. */
 if(rho<1e-10*s.b){
  const double eps=std::max(1e-4*(s.local.hole[0].mass+s.local.hole[1].mass),1e-8*s.b);
  for(int v=0;v<4;v++)for(int d=0;d<4;d++)f[v][d]=0;
  for(int k=0;k<4;k++){
   double xp[3]={x[0],0,0};xp[1+k/2]=(k%2? -eps:eps);Fields fp{};sample_fields(s,xp,fp);
   for(int v=0;v<4;v++)for(int d=0;d<4;d++)f[v][d]+=.25*fp[v][d];
  }
  /* Values have an exact map-axis limit: the Fourier zero mode at the
   * central A,B coordinates. Averaging displaced values introduces an
   * O(eps^2) bias that an independent Cartesian Hessian amplifies as h^-2.
   * Retain the transverse approximation only for Cartesian derivatives. */
  double A=-1,B;
  if(std::abs(x[0])<=s.b){
   const double R=std::acos(std::clamp(x[0]/s.b,-1.0,1.0));
   B=std::tan(R/2-Piq);
  }else{
   const double X=std::acosh(std::abs(x[0])/s.b);
   A=2*std::tanh(X/2)-1;B=x[0]>0?-1:1;
  }
  std::vector<double>Ta,Da,Tb,Db;
  chebyshev_basis(s.local.n[0],A,Ta,Da);chebyshev_basis(s.local.n[1],B,Tb,Db);
  for(int v=0;v<4;v++){
   f[v][0]=0;
   for(int j=0;j<s.local.n[1];j++)for(int i=0;i<s.local.n[0];i++)
    f[v][0]+=.5*(A-1)*s.coefficients[v+4*(i+s.local.n[0]*j)]*Ta[i]*Tb[j];
  }
  return;
 }
 /* Stable inverse of x=b cosh(X) cos(R), rho=b sinh(X) sin(R).
  * Rationalize the small root rather than subtracting nearly equal
  * puncture distances (which loses R close to the axis). */
 const double w0=x[0]/s.b,q=rho/s.b,t=.5*(w0*w0+q*q-1),root=std::hypot(t,q);
 const double sx2=t<0?q*q/(root-t):root+t;
 const double sr2=t>0?q*q/(root+t):root-t;
 double X=std::asinh(std::sqrt(sx2)),R=std::atan2(std::sqrt(sr2),w0/std::sqrt(1+sx2));
 double A=2*std::tanh(X/2)-1,B=std::tan(R/2-Piq),phi=std::atan2(x[2],x[1]);if(phi<0)phi+=2*Pi;
 double dat[10][4]={};derivs w{};double **ptr[10]={&w.d0,&w.d1,&w.d2,&w.d3,&w.d11,&w.d12,&w.d13,&w.d22,&w.d23,&w.d33};
 for(int d=0;d<10;d++)*ptr[d]=dat[d];
 /* Differentiate the coefficient basis itself, including the cosine
  * Nyquist derivative which vanishes only at collocation phi values. */
 double V[4][4];evaluate_coefficients(s,A,B,phi,V);
 for(int v=0;v<4;v++){
  dat[0][v]=(A-1)*V[v][0];dat[1][v]=V[v][0]+(A-1)*V[v][1];dat[2][v]=(A-1)*V[v][2];dat[3][v]=(A-1)*V[v][3];
 }
 double xx,rr,y,z;AB_To_XR(4,A,B,&X,&R,&w);C_To_c(4,X,R,&xx,&rr,s.b,&w);
 rx3_To_xyz(4,xx,rr,phi,&y,&z,&w);for(int v=0;v<4;v++)for(int d=0;d<4;d++)f[v][d]=dat[d][v];
}
}
extern "C" {
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
  if(sampler_only)return s;
  s->geometry.resize(s->npt);allocate_derivs(&s->work,s->ntotal);
  s->derivatives.initialize(c->n);
  for(int k=0;k<c->n[2];k++)for(int j=0;j<c->n[1];j++)for(int i=0;i<c->n[0];i++){
   Fields f{};double x[3];transform(*s,i,j,k,f,x);hispid::Background bg;hispid::background(s->local,x,bg);
   int p=pindex(*s,i,j,k);cache(bg,s->geometry[p]);
   double sn=std::sin(Pih*(2*i+1)/c->n[0])*std::sin(Pih*(2*j+1)/c->n[1]);s->geometry[p].weight=std::pow(sn,6);
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
 try{
  residual(*s,s->values.data(),r.data());
  for(int it=0;it<=s->local.max_newton;it++){
   update_diag(*s,r);double err=norminf(r);if(err<=s->local.tolerance){s->diag.converged=1;break;}
   if(it==s->local.max_newton)break;
   s->diag.newton_iterations++;
   Sparse M=preconditioner(*s);for(int i=0;i<s->ntotal;i++)rhs[i]=-r[i];
   double forcing=std::min(.05,std::max(1e-5,std::sqrt(err)));
   if(!gmres(*s,M,rhs,step,forcing)){hispid::last_error="Krylov iteration limit";break;}
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
static int sample_with_derivatives(HiSpID_Data*s,int count,const double*xyz,HiSpID_Point*out,double*dgamma){
 if(!s||!xyz||!out||count<0)return -1;
 try{for(int p=0;p<count;p++){
  double x[3];to_local(*s,xyz+3*p,x);hispid::Background bg;hispid::background(s->local,x,bg);Fields f{};sample_fields(*s,x,f);
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
 std::fill(out,out+7,0.0);double delta=radius*1e-4;
 /* Gauss-Legendre cos(theta) and uniform phi. */
 for(int a=0;a<nt;a++){
  double mu=std::cos(Pi*(a+.75)/(nt+.5)),pp=0;
  for(int it=0;it<20;it++){
   double p0=1,p1=mu;for(int l=2;l<=nt;l++){double p=((2*l-1)*mu*p1-(l-1)*p0)/l;p0=p1;p1=p;}
   pp=nt*(mu*p1-p0)/(mu*mu-1);double step=p1/pp;mu-=step;if(std::abs(step)<1e-15)break;
  }double weight=2/((1-mu*mu)*pp*pp)*(2*Pi/np)*radius*radius;
  for(int b=0;b<np;b++){
   double ph=2*Pi*(b+.5)/np,n[3]={std::sqrt(1-mu*mu)*std::cos(ph),std::sqrt(1-mu*mu)*std::sin(ph),mu};
   double x[3];for(int i=0;i<3;i++)x[i]=center[i]+radius*n[i];HiSpID_Point p;
   if(HiSpID_sample(s,1,x,&p))return -2;
   double dg[3][9];for(int d=0;d<3;d++){
    double xx[12];HiSpID_Point st[4];int offsets[4]={-2,-1,1,2};for(int k=0;k<4;k++)for(int i=0;i<3;i++)xx[3*k+i]=x[i]+(i==d?offsets[k]*delta:0);
    if(HiSpID_sample(s,4,xx,st))return -2;
    for(int ij=0;ij<9;ij++)dg[d][ij]=(st[0].gamma[ij]-8*st[1].gamma[ij]+8*st[2].gamma[ij]-st[3].gamma[ij])/(12*delta);
   }
   double E=0,P[3]={};for(int i=0;i<3;i++)for(int j=0;j<3;j++){
    E+=n[i]*(dg[j][3*i+j]-dg[i][3*j+j]);
    P[i]+=(p.Kij[3*i+j]-p.mean_curvature*p.gamma[3*i+j])*n[j];
   }
   out[0]+=weight*E/(16*Pi);for(int i=0;i<3;i++){
    out[i+1]+=weight*P[i]/(8*Pi);out[i+4]+=weight*radius*(n[(i+1)%3]*P[(i+2)%3]-n[(i+2)%3]*P[(i+1)%3])/(8*Pi);
   }
  }
 }return 0;
}
void HiSpID_destroy(HiSpID_Data*s){if(!s)return;if(s->work)free_derivs(s->work);delete s;}
}
