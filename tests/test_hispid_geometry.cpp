#include "../src/HiSpID_internal.hpp"
#include "../src/HiSpID_spectral.hpp"
#include <cstdio>
#include <cstdlib>
using namespace hispid;
static void close(double a,double b,double tol,const char*label){
 if(std::abs(a-b)>tol){std::fprintf(stderr,"%s: %.17g vs %.17g\n",label,a,b);std::exit(1);}
}
int main(){
 const int shape[3]={12,10,8},size=4*shape[0]*shape[1]*shape[2];
 derivs *legacy,*matrix;allocate_derivs(&legacy,size);allocate_derivs(&matrix,size);
 unsigned state=12345;
 for(int p=0;p<size;p++){state=1664525*state+1013904223;legacy->d0[p]=matrix->d0[p]=(double)state/4294967296.0-.5;}
 Derivatives_AB3(4,shape[0],shape[1],shape[2],legacy);
 SpectralDerivatives spectral;spectral.initialize(shape);spectral.apply(4,matrix);
 double *ld[9]={legacy->d1,legacy->d2,legacy->d3,legacy->d11,legacy->d12,legacy->d13,legacy->d22,legacy->d23,legacy->d33};
 double *md[9]={matrix->d1,matrix->d2,matrix->d3,matrix->d11,matrix->d12,matrix->d13,matrix->d22,matrix->d23,matrix->d33};
 for(int d=0;d<9;d++){
  double error=0,scale=1;for(int p=0;p<size;p++){error=std::max(error,std::abs(ld[d][p]-md[d][p]));scale=std::max(scale,std::abs(ld[d][p]));}
  close(error/scale,0,1e-11,"cached differentiation matrices versus legacy interpolant");
 }
 free_derivs(legacy);free_derivs(matrix);
 Jet x=Jet::variable(.7,1),y=Jet::variable(-.4,2),z=Jet::variable(.3,3);
 Jet w=Jet(.03)*x*x+Jet(.02)*y+Jet(.01)*z*z;
 Jet g[3][3],inv[3][3],C[3][3][3];
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)g[i][j]=i==j?exp(Jet(4)*w):Jet(0);
 invert(g,inv);connection(g,inv,C);
 long double grad2=0,lapw=0;for(int i=1;i<4;i++){grad2+=w.d[i]*w.d[i];lapw+=w.h[i][i];}
 close(curvature(inv,C),std::exp(-4*w.v)*(-8*lapw-8*grad2),1e-12,"manufactured curved Ricci");
 Jet f=x*x+Jet(2)*y*y+Jet(3)*z*z;long double grad=0;for(int i=1;i<4;i++)grad+=w.d[i]*f.d[i];
 close(laplacian(inv,C,f),std::exp(-4*w.v)*(12+2*grad),1e-12,"manufactured curved scalar Laplacian");
 Jet bv[3]={1,2,3},L[3][3];double div[3];
 longitudinal(g,C,bv,L);divergence(inv,C,L,div);for(int i=0;i<3;i++)close(div[i],0,1e-12,"curved translation Killing");
 bv[0]=x;bv[1]=y;bv[2]=z;longitudinal(g,C,bv,L);divergence(inv,C,L,div);for(int i=0;i<3;i++)close(div[i],0,1e-12,"curved dilation Killing");
 bv[0]=-y;bv[1]=x;bv[2]=0;longitudinal(g,C,bv,L);divergence(inv,C,L,div);for(int i=0;i<3;i++)close(div[i],0,1e-12,"curved rotation Killing");
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)g[i][j]=Jet(i==j?1:0);
 invert(g,inv);connection(g,inv,C);bv[0]=x*x;bv[1]=bv[2]=0;
 longitudinal(g,C,bv,L);divergence(inv,C,L,div);close(div[0],8.0/3,1e-12,"flat longitudinal vector Laplacian");
 HiSpID_Hole hole{};hole.mass=1;hole.spin[0]=.2;hole.spin[1]=.3;hole.spin[2]=.4;
 hole.velocity[0]=.2;hole.velocity[1]=-.1;hole.velocity[2]=.15;
 double point[3]={1.1,-.7,.4};Seed s;seed(hole,1,point,s);double h=1e-3;
 for(int d=0;d<3;d++){
  Seed st[4];int offset[4]={-2,-1,1,2};for(int k=0;k<4;k++){double p[3]={point[0],point[1],point[2]};p[d]+=offset[k]*h;seed(hole,1,p,st[k]);}
  for(int i=0;i<3;i++)for(int j=0;j<3;j++){
   long double dg=(st[0].physical[i][j].v-8*st[1].physical[i][j].v+8*st[2].physical[i][j].v-st[3].physical[i][j].v)/(12*h);
   long double dK=(st[0].extrinsic[i][j].v-8*st[1].extrinsic[i][j].v+8*st[2].extrinsic[i][j].v-st[3].extrinsic[i][j].v)/(12*h);
   long double ddg=(-st[0].physical[i][j].v+16*st[1].physical[i][j].v-30*s.physical[i][j].v+16*st[2].physical[i][j].v-st[3].physical[i][j].v)/(12*h*h);
   close(dg,s.physical[i][j].d[d+1],1e-10,"seed metric AD derivative");
   close(dK,s.extrinsic[i][j].d[d+1],1e-10,"seed extrinsic AD derivative");
   close(ddg,s.physical[i][j].h[d+1][d+1],1e-8,"seed metric AD Hessian");
  }
 }
 std::puts("manufactured curvature/scalar/vector operators and independent seed derivatives passed");return 0;
}
