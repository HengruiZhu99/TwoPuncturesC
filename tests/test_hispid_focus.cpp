#include "../src/HiSpID_modal_block.hpp"
#include <cstdio>
int main(){
 gsl_set_error_handler_off();double worst=0,compat=0;
 for(int m=0;m<=4;m++){
  const int a=96,b=144,N=a*b;auto w=hispid::focus_weights(a,b,m);
  for(int side=0;side<hispid::focus_rank;side++){
   long double sum=0,physical=0;for(int j=0;j<b;j++)for(int i=0;i<a;i++){
    double z=-std::cos(std::acos(-1.)*(i+.5)/a),s=(1+z)/2,t=hispid::AxisDerivatives::radial_t(s),eta=hispid::AxisDerivatives::angular_eta(-std::cos(std::acos(-1.)*(j+.5)/b));
    // Independent Cartesian-smoothness polynomial at both foci.
    double P=1+(m+1)*t+.5*(m+1)*(m+2)*t*t+.2*eta*(1+(m+3)*t+.5*(m+2)*(m+5)*t*t);sum+=(long double)w[side*N+i+a*j]*P;
    // Independently prescribed smooth Cartesian field:
    // Re[(y+i*z)^m]/(1+x*x+y*y+z*z)^(m+2), at separation half-width .25.
    const double q=1-t,den=q*q*(1+.0625*eta*eta)+.25*t;
    const double cartesian_P=-.5*std::pow(.5,m)*std::pow(q,m+3)/std::pow(den,m+2);
    physical+=(long double)w[side*N+i+a*j]*cartesian_P;
   }compat=std::max({compat,double(std::abs(sum)),double(std::abs(physical))});
  }
  const int na=6,nb=9,n=na*nb;hispid::ModalBlock B;B.na=na;B.nb=nb;B.focus_endpoint=hispid::focus_weights(na,nb,m);
  B.lu.assign(na*na*nb,0);B.transfer.resize(na*na*nb);B.lower.assign(n,0);B.upper.assign(n,0);B.permutation.resize(n);
  std::vector<double>A(n*n,0),f(n),x(4*n,0);
  for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
   int r=i+na*j;A[r*n+r]=B.lu[(j*na+i)*na+i]=5+.1*i;
   for(int d:{-1,1})if(i+d>=0&&i+d<na)A[r*n+r+d]=B.lu[(j*na+i)*na+i+d]=-.3;
   if(j)A[r*n+r-na]=B.lower[r]=-.2;if(j+1<nb)A[r*n+r+na]=B.upper[r]=-.4;
   f[r]=std::sin(.3*r);x[4*r]=f[r];
  }
  for(int side=0;side<hispid::focus_rank;side++){int r=hispid::focus_row(side,na,nb);for(int p=0;p<n;p++)A[r*n+p]=B.focus_endpoint[side*n+p];}
  B.factor();B.solve(x.data(),0,0);
  auto mat=gsl_matrix_view_array(A.data(),n,n);auto rhs=gsl_vector_view_array(f.data(),n);auto perm=gsl_permutation_alloc(n);int sign;
  if(gsl_linalg_LU_decomp(&mat.matrix,perm,&sign)||gsl_linalg_LU_svx(&mat.matrix,perm,&rhs.vector))return 1;
  gsl_permutation_free(perm);for(int p=0;p<n;p++)worst=std::max(worst,std::abs(x[4*p]-f[p])/(1+std::abs(f[p])));
 }
 printf("focus polynomial compatibility %.17g; independent dense inverse %.17g\n",compat,worst);
 return compat<1e-12&&worst<1e-10?0:1;
}
