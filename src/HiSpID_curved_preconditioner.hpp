#ifndef HISPID_CURVED_PRECONDITIONER_HPP
#define HISPID_CURVED_PRECONDITIONER_HPP
#include "HiSpID_axis.hpp"
#include "HiSpID_cache_kernels.hpp"
namespace hispid {
// Curved scalar-Laplacian symbol with mixed meridional derivatives.
// The vector m=1 symbol projects the complete operator onto axial rotation;
// other vector modes use the mean principal part,10/9 times Laplace. Full
// component couplings always remain in the fine-grid vector block.
struct ModalSymbol { double zero=0,first[2]={},second[2]={},mixed=0; };
HISPID_GEOMETRY_INLINE ModalSymbol curved_symbol(const Cached&g,double b,
        double axial,double y,double z,int mode,int group){
 using J=JetT<double>;
 J x=J::variable(axial,1),Y=J::variable(y,2),Z=J::variable(z,3);
 J rho2=Y*Y+Z*Z,plus=sqrt((x+J(b))*(x+J(b))+rho2),minus=sqrt((x-J(b))*(x-J(b))+rho2);
 J C=(plus+minus)/J(2*b),eta=(plus-minus)/J(2*b),t=(C-J(1))/(C+J(1));
 J sigma=AxisDerivatives::radial_inverse(t);
 const double k=AxisDerivatives::angular_stretch,T=std::tanh(k);
 J e=eta*J(T),zeta=unary(e,std::atanh(e.v),1/(1-e.v*e.v),2*e.v/std::pow(1-e.v*e.v,2))/J(k);
 int r=AxisDerivatives::exponent(mode);
 // Use Cartesian harmonic polynomials for the regular r<=4 factor.
 // Differentiating rho^r and exp(i*m*phi) separately cancels O(rho^-2)
 // terms for m=r and loses precision near the axis.
 J radial=J(-2)*(J(1)-t),real(1),imag;
 if(r)radial=radial/power(J(b)*(C+J(1)),r);
 for(int n=0;n<r;n++){J next=real*Y-imag*Z;imag=real*Z+imag*Y;real=next;}
 real=radial*real;imag=radial*imag;
 const double angle=r*std::atan2(z,y),co=std::cos(angle),si=std::sin(angle);
 const double denominator=radial.v*std::pow(std::sqrt(rho2.v),r);
 const double pg[3]={0,-z/rho2.v,y/rho2.v},p2=rho2.v*rho2.v;
 const double ph[3][3]={{0,0,0},{0,2*y*z/p2,(z*z-y*y)/p2},{0,(z*z-y*y)/p2,-2*y*z/p2}};
 const int remaining=mode-r;
 auto apply=[&](const J&fr,const J&fi){
  double q=0;
  for(int d=0;d<10;d++){
   double vr=fr.v,vi=fi.v;
   if(d>0&&d<4){vr=fr.d[d]-remaining*pg[d-1]*fi.v;vi=fi.d[d]+remaining*pg[d-1]*fr.v;}
   if(d>=4){
    int a=hessian_i(d-4),b=hessian_j(d-4);
    vr=fr.h[a][b]-remaining*(pg[a-1]*fi.d[b]+pg[b-1]*fi.d[a]+ph[a-1][b-1]*fi.v)-remaining*remaining*pg[a-1]*pg[b-1]*fr.v;
    vi=fi.h[a][b]+remaining*(pg[a-1]*fr.d[b]+pg[b-1]*fr.d[a]+ph[a-1][b-1]*fr.v)-remaining*remaining*pg[a-1]*pg[b-1]*fi.v;
   }
   if(group&&mode==1){
    // beta=(0,-Im(f),Re(f)): its Cartesian m=1 pair is the axisymmetric
    // toroidal field. Keep the cross-component connections/curvature
    // responsible for the axial Killing identity in this difficult sector.
    const double My=-g.vec[1][1][d]*vi+g.vec[1][2][d]*vr;
    const double Mz=-g.vec[2][1][d]*vi+g.vec[2][2][d]*vr;
    q+=-si*My+co*Mz;
   }else q+=(group?10./9:1.)*g.lap[d]*(co*vr+si*vi);
  }
  return q/denominator;
 };
 ModalSymbol out;out.zero=apply(real,imag);
 const J coordinates[2]={sigma-J(sigma.v),zeta-J(zeta.v)};
 for(int a=0;a<2;a++){
  out.first[a]=apply(real*coordinates[a],imag*coordinates[a]);
  out.second[a]=apply(real*coordinates[a]*coordinates[a]/J(2),imag*coordinates[a]*coordinates[a]/J(2));
 }
 out.mixed=apply(real*coordinates[0]*coordinates[1],imag*coordinates[0]*coordinates[1]);
 return out;
}
}
#endif
