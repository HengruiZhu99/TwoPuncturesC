#ifndef HISPID_JETS_HPP
#define HISPID_JETS_HPP
#include <cmath>
#include <stdexcept>
namespace hispid {
/* Second-order forward jets. Long double protects the signed-lapse
 * cancellations near the QI throat. No finite differences in seed geometry. */
struct Jet {
  long double v, d[4], h[4][4];
  Jet(long double x=0):v(x),d{},h{} {}
  static Jet variable(long double v,int k){Jet a(v);a.d[k]=1;return a;}
};
inline Jet operator+(const Jet&a,const Jet&b){
  Jet c(a.v+b.v);for(int i=0;i<4;i++){c.d[i]=a.d[i]+b.d[i];
    for(int j=0;j<4;j++)c.h[i][j]=a.h[i][j]+b.h[i][j];}return c;
}
inline Jet operator-(const Jet&a,const Jet&b){
  Jet c(a.v-b.v);for(int i=0;i<4;i++){c.d[i]=a.d[i]-b.d[i];
    for(int j=0;j<4;j++)c.h[i][j]=a.h[i][j]-b.h[i][j];}return c;
}
inline Jet operator-(const Jet&a){return Jet(0)-a;}
inline Jet operator*(const Jet&a,const Jet&b){
  Jet c(a.v*b.v);for(int i=0;i<4;i++){c.d[i]=a.d[i]*b.v+a.v*b.d[i];
    for(int j=0;j<4;j++)c.h[i][j]=a.h[i][j]*b.v+a.d[i]*b.d[j]
      +a.d[j]*b.d[i]+a.v*b.h[i][j];}return c;
}
inline Jet unary(const Jet&a,long double v,long double p,long double pp){
  Jet c(v);for(int i=0;i<4;i++){c.d[i]=p*a.d[i];
    for(int j=0;j<4;j++)c.h[i][j]=p*a.h[i][j]+pp*a.d[i]*a.d[j];}return c;
}
inline Jet power(const Jet&a,long double p){
  long double v=std::pow(a.v,p);return unary(a,v,p*v/a.v,p*(p-1)*v/(a.v*a.v));
}
inline Jet inverse(const Jet&a){return power(a,-1);}
inline Jet operator/(const Jet&a,const Jet&b){return a*inverse(b);}
inline Jet sqrt(const Jet&a){return power(a,.5L);}
inline Jet exp(const Jet&a){long double v=std::exp(a.v);return unary(a,v,v,v);}
inline Jet tan(const Jet&a){long double v=std::tan(a.v),p=1+v*v;return unary(a,v,p,2*v*p);}
inline Jet tanh(const Jet&a){long double v=std::tanh(a.v),p=1-v*v;return unary(a,v,p,-2*v*p);}
/* This differentiated jet is correct to FIRST order only. Used for K and
 * connection first derivatives; third spacetime derivatives are never used. */
inline Jet diff(const Jet&a,int k){Jet c(a.d[k]);for(int i=0;i<4;i++)c.d[i]=a.h[k][i];return c;}
inline Jet determinant(const Jet g[3][3]){
 return g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[2][1])
       -g[0][1]*(g[1][0]*g[2][2]-g[1][2]*g[2][0])
       +g[0][2]*(g[1][0]*g[2][1]-g[1][1]*g[2][0]);
}
inline void invert(const Jet g[3][3],Jet inv[3][3]){
 Jet det=determinant(g);if(!(det.v>0))throw std::runtime_error("nonpositive spatial determinant");
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  int a=(j+1)%3,b=(j+2)%3,c=(i+1)%3,d=(i+2)%3;
  inv[i][j]=(g[a][c]*g[b][d]-g[a][d]*g[b][c])/det;
 }
}
inline void connection(const Jet g[3][3],const Jet inv[3][3],Jet C[3][3][3]){
 for(int k=0;k<3;k++)for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  C[k][i][j]=0;for(int l=0;l<3;l++)C[k][i][j]=C[k][i][j]
   +Jet(.5L)*inv[k][l]*(diff(g[j][l],i+1)+diff(g[i][l],j+1)-diff(g[i][j],l+1));
 }
}
inline long double curvature(const Jet inv[3][3],const Jet C[3][3][3]){
 long double R=0;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){
  long double rij=0;for(int k=0;k<3;k++){
   rij+=C[k][i][j].d[k+1]-C[k][i][k].d[j+1];
   for(int l=0;l<3;l++)rij+=C[k][i][j].v*C[l][k][l].v-C[l][i][k].v*C[k][j][l].v;
  }R+=inv[i][j].v*rij;
 }
 return R;
}
}
#endif
