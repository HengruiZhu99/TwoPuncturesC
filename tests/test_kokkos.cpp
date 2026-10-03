/* Manufactured nonsymmetric solve and independent polynomial/coordinate
 * controls, executed unchanged on Serial, OpenMP and CUDA. */
#include "PunctureKokkos.hpp"
#include "HiSpID_spectral.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
using namespace puncture;
static int checks;
#define CHECK(c) do{checks++;if(!(c)){fprintf(stderr,"line%d: %s\n",__LINE__,#c);return 1;}}while(0)
int main(){
 initialize();{
 const int n=19;std::vector<double>matrix(n*n,0),exact(n),rhs(n,0),answer(n);
 for(int i=0;i<n;i++){exact[i]=std::sin(i+.2);for(int j=0;j<n;j++)matrix[n*i+j]=i==j?4+.1*i:(j==i+1?-.8:(i==j+2?.35:0));}
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)rhs[i]+=matrix[n*i+j]*exact[j];
 auto a=upload(matrix,"matrix"),b=upload(rhs,"RHS");View x("solution",n);
 Action A=[=](View in,View out){Kokkos::parallel_for("dense A",Range(0,n),KOKKOS_LAMBDA(int i){double sum=0;for(int j=0;j<n;j++)sum+=a(n*i+j)*in(j);out(i)=sum;});};
 Action M=[=](View in,View out){each(out,KOKKOS_LAMBDA(int i){out(i)=in(i)/a(n*i+i);});};
 for(int method=0;method<2;method++)for(int restart:{3,19}){
  PK_Options o{method,200,restart,1,0,0,1e-11,nullptr,nullptr};PK_Result r{};zero(x);
  CHECK(solve(b,x,o,A,M,r)==PK_SUCCESS);CHECK(r.true_residual<=o.absolute_tolerance);download(x,answer.data());
  double residual=0;for(int i=0;i<n;i++){double d=rhs[i];for(int j=0;j<n;j++)d-=matrix[n*i+j]*answer[j];residual+=d*d;CHECK(std::abs(answer[i]-exact[i])<1e-10);}CHECK(std::sqrt(residual)<=o.absolute_tolerance);
  o.max_iterations=1;o.absolute_tolerance=1e-15;zero(x);CHECK(solve(b,x,o,A,M,r)==PK_LIMIT);CHECK(r.iterations==1);
  Action fail=[](View,View){throw std::runtime_error("manufactured failure");};zero(x);CHECK(solve(b,x,o,A,fail,r)==PK_CALLBACK);
  Action nan=[](View,View out){Kokkos::deep_copy(out,NAN);};zero(x);CHECK(solve(b,x,o,A,nan,r)==PK_NONFINITE);
 }
 PK_Options valid{PK_GMRES,200,19,1,0,0,1e-11,nullptr,nullptr};PK_Result r{};
 Action identity=[](View in,View out){copy(out,in);};Action empty;
 View z("zero RHS",n);zero(z);
 for(int method:{PK_GMRES,PK_BICGSTAB}){
  valid.method=method;
  for(auto rhs_view:{b,z}){zero(x);CHECK(solve(rhs_view,x,valid,empty,M,r)==PK_INVALID);CHECK(solve(rhs_view,x,valid,A,empty,r)==PK_INVALID);}
  zero(x);CHECK(solve(z,x,valid,identity,identity,r)==PK_SUCCESS);CHECK(r.iterations==0&&r.true_residual==0);
  copy(x,b);CHECK(solve(b,x,valid,identity,identity,r)==PK_SUCCESS);CHECK(r.true_residual<=valid.absolute_tolerance);
  if(method==PK_BICGSTAB)CHECK(r.iterations==0); // nonzero initial guess is retained
  for(int bad=0;bad<8;bad++){
   auto o=valid;switch(bad){case 0:o.method=2;break;case 1:o.max_iterations=0;break;case 2:o.restart=4097;break;case 3:o.verify_true_residual=2;break;case 4:o.legacy_bicgstab=2;break;case 5:o.eager_gmres_basis=2;break;case 6:o.absolute_tolerance=NAN;break;case 7:o.absolute_tolerance=-1;}
   zero(x);CHECK(solve(b,x,o,A,M,r)==PK_INVALID);
  }
  zero(x);CHECK(solve(b,b,valid,A,M,r)==PK_INVALID);
 }
 valid.method=PK_BICGSTAB;valid.verify_true_residual=0;zero(x);CHECK(solve(b,x,valid,A,M,r)==PK_INVALID);
 valid.verify_true_residual=1;Action singular=[](View,View out){zero(out);};CHECK(solve(b,x,valid,singular,identity,r)==PK_BREAKDOWN);
 int shape[3]={12,16,10};Operator op(shape,1,false);initialize_spectral(op,shape,2.);
 std::vector<double>v(12*16*10,1),out(v.size());auto in=upload(v,"constant");op.spectral.apply(in);
 for(int d=1;d<10;d++){download(op.spectral.work[d],out.data());for(double z:out)CHECK(z==0);}
 // Polynomial and Fourier jets, including the even-grid Nyquist second
 // derivative, are compared with analytic values on the actual device.
 for(int nyquist=0;nyquist<2;nyquist++){
  std::vector<double>jets(v.size()*10);
  for(int k=0;k<shape[2];k++)for(int j=0;j<shape[1];j++)for(int i=0;i<shape[0];i++){
   int p=i+shape[0]*(j+shape[1]*k);double a=-std::cos(Pi*(i+.5)/shape[0]),b0=-std::cos(Pi*(j+.5)/shape[1]),phi=2*Pi*k/shape[2];auto*f=jets.data()+10*p;
   if(nyquist){f[0]=k%2?-1:1;f[9]=-25*f[0];}
   else{double co=std::cos(2*phi),si=std::sin(2*phi);f[0]=a*a+3*b0*b0*b0+(a+b0)*co;f[1]=2*a+co;f[2]=9*b0*b0+co;f[3]=-2*(a+b0)*si;f[4]=2;f[6]=f[8]=-2*si;f[7]=18*b0;f[9]=-4*(a+b0)*co;}
   v[p]=f[0];
  }
  in=upload(v,"analytic spectral control");op.spectral.apply(in);
  for(int d=0;d<10;d++){download(op.spectral.work[d],out.data());for(size_t p=0;p<v.size();p++)CHECK(std::abs(out[p]-jets[10*p+d])/(1+std::abs(jets[10*p+d]))<2e-11);}
  // Isolate the device coordinate chain from spectral rounding near the
  // singular prolate foci: feed the same computed jets to the independent
  // GSL chain. The analytic derivative gate above remains separate.
  for(int d=0;d<10;d++){download(op.spectral.work[d],out.data());for(size_t p=0;p<v.size();p++)jets[10*p+d]=out[p];}
  op.physical_fields(in);std::vector<double>cartesian(jets.size());download(op.fields,cartesian.data());
  for(int k=0;k<shape[2];k++)for(int j=0;j<shape[1];j++)for(int i=0;i<shape[0];i++){
   int p=i+shape[0]*(j+shape[1]*k);double a=-std::cos(Pi*(i+.5)/shape[0]),b0=-std::cos(Pi*(j+.5)/shape[1]),phi=2*Pi*k/shape[2],f[10];std::memcpy(f,jets.data()+10*p,sizeof(f));double original[10];std::memcpy(original,f,sizeof(f));for(int d=0;d<10;d++)f[d]*=a-1;f[1]+=original[0];f[4]+=2*original[1];f[5]+=original[2];f[6]+=original[3];
   derivs w{};double**ptr[10]={&w.d0,&w.d1,&w.d2,&w.d3,&w.d11,&w.d12,&w.d13,&w.d22,&w.d23,&w.d33};for(int d=0;d<10;d++)*ptr[d]=f+d;
   double X,R,x0,r0,y,z0;AB_To_XR(1,a,b0,&X,&R,&w);C_To_c(1,X,R,&x0,&r0,2.,&w);rx3_To_xyz(1,x0,r0,phi,&y,&z0,&w);
   for(int d=0;d<10;d++)CHECK(std::abs(cartesian[10*p+d]-f[d])/(1+std::abs(f[d]))<1e-10);
  }
 }
 // Coordinate chain compared against the inherited GSL implementation on
 // every collocation point, with all ten independent input jets.
 auto chain=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),op.chain);
 for(int k=0;k<shape[2];k++)for(int j=0;j<shape[1];j++)for(int i=0;i<shape[0];i++){
  double A0=-std::cos(Pi*(i+.5)/shape[0]),B=-std::cos(Pi*(j+.5)/shape[1]),phi=2*Pi*k/shape[2];
  for(int basis=0;basis<10;basis++){
   double values[10]{},portable[10]{};values[basis]=portable[basis]=1;
   coordinate_chain(A0,B,chain.data()+6*(i+shape[0]*j),std::cos(phi),std::sin(phi),portable);
   double original[10];std::memcpy(original,values,sizeof(values));for(int d=0;d<10;d++)values[d]*=A0-1;values[1]+=original[0];values[4]+=2*original[1];values[5]+=original[2];values[6]+=original[3];
   derivs w{};double**ptr[10]={&w.d0,&w.d1,&w.d2,&w.d3,&w.d11,&w.d12,&w.d13,&w.d22,&w.d23,&w.d33};for(int d=0;d<10;d++)*ptr[d]=values+d;
   double X,R,x,r,y,z;AB_To_XR(1,A0,B,&X,&R,&w);C_To_c(1,X,R,&x,&r,2.,&w);rx3_To_xyz(1,x,r,phi,&y,&z,&w);
   for(int d=0;d<10;d++)CHECK(std::abs(values[d]-portable[d])/(1+std::abs(values[d]))<3e-13);
  }
 }
 }printf("Kokkos %s independent controls: %d passed\n",Exec::name(),checks);return 0;
}
