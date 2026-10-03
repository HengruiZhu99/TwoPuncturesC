#include "PunctureKokkos.hpp"
#include "HiSpID_axis.hpp"
#include <utility>
namespace puncture {
void initialize_hispid_spectral(Operator&o,hispid::AxisDerivatives&axis,const int*n,double b){
 std::copy(n,n+3,axis.n.begin());axis.raw.n=axis.n;
 View first[2],second[2];
 for(int a=0;a<3;a++){
  const int N=n[a];auto D=View("D kernel cache",N*N),D2=View("D2 kernel cache",N*N);
  o.spectral.D[a]=D;o.spectral.D2[a]=D2;
  if(a<2){
   View nodes("barycentric nodes",2*N);
   Kokkos::parallel_for("barycentric nodes",Range(0,N),KOKKOS_LAMBDA(int i){double t=Pi*(i+.5)/N;nodes(2*i)=-std::cos(t);nodes(2*i+1)=(i%2?-1:1)*std::sin(t);});
   Kokkos::parallel_for("first derivative cache",Range(0,N),KOKKOS_LAMBDA(int i){double diagonal=0;
    for(int j=0;j<N;j++)if(j!=i){double value=nodes(2*j+1)/(nodes(2*i+1)*(nodes(2*i)-nodes(2*j)));D(i*N+j)=value;diagonal-=value;}D(i*N+i)=diagonal;
   });
   Kokkos::parallel_for("second derivative cache",Range(0,N*N),KOKKOS_LAMBDA(int q){int i=q/N,j=q%N;double value=0;for(int k=0;k<N;k++)value+=D(i*N+k)*D(k*N+j);D2(q)=value;});
   Kokkos::parallel_for("constant derivative identity",Range(0,N),KOKKOS_LAMBDA(int i){double diagonal=0;for(int j=0;j<N;j++)if(j!=i)diagonal-=D2(i*N+j);D2(i*N+i)=diagonal;});
   auto coordinate=View("mapped nodes",N),ma=View("first map",N),mma=View("second map",N),coefficient=View("Chebyshev coefficient cache",N*N);const int direction=a;
   o.spectral.coordinate[a]=coordinate;o.spectral.coefficient[a]=coefficient;first[a]=ma;second[a]=mma;
   Kokkos::parallel_for("mapped coordinate cache",Range(0,N),KOKKOS_LAMBDA(int i){double values[3];hispid::AxisDerivatives::node(direction,N,i,values);coordinate(i)=values[0];ma(i)=values[1];mma(i)=values[2];});
   Kokkos::parallel_for("Chebyshev coefficient cache",Range(0,N*N),KOKKOS_LAMBDA(int q){coefficient(q)=hispid::AxisDerivatives::coefficient_at(N,q/N,q%N);});
   for(auto pair:{std::make_pair(coordinate,&axis.coordinate[a]),std::make_pair(ma,&axis.first_map[a]),std::make_pair(mma,&axis.second_map[a]),std::make_pair(coefficient,&axis.coefficient[a])}){
    pair.second->resize(pair.first.extent(0));download(pair.first,pair.second->data());
   }
  }else{
   Kokkos::parallel_for("Fourier derivative caches",Range(0,N*N),KOKKOS_LAMBDA(int q){int i=q/N,j=q%N;
    if(i==j){D(q)=0;D2(q)=-(N*N+2.0)/12;}else{double t=Pi*(i-j)/N,sign=(std::abs(i-j)%2?-1:1);D(q)=.5*sign/std::tan(t);D2(q)=-.5*sign/(std::sin(t)*std::sin(t));}
   });
  }
  axis.raw.D[a].resize(N*N);axis.raw.D2[a].resize(N*N);download(D,axis.raw.D[a].data());download(D2,axis.raw.D2[a].data());
 }
 const int na=n[0],nb=n[1],np=n[2],stride=na*nb,points=stride*np;
 auto forward=View("Fourier forward cache",np*np),inverse=View("Fourier inverse cache",np*np),phi=View("Fourier first cache",np*np),phi2=View("Fourier second cache",np*np);
 o.spectral.inverse=inverse;o.spectral.inverse_phi=phi;o.spectral.inverse_phi2=phi2;
 Kokkos::parallel_for("Fourier partner caches",Range(0,np*np),KOKKOS_LAMBDA(int q){int k=q/np,mode=q%np;double values[3];hispid::AxisDerivatives::fourier(np,k,mode,values);forward(mode*np+k)=values[0];inverse(q)=values[0];phi(q)=values[1];phi2(q)=values[2];});
 for(auto pair:{std::make_pair(forward,&axis.forward),std::make_pair(inverse,&axis.inverse),std::make_pair(phi,&axis.inverse_phi),std::make_pair(phi2,&axis.inverse_phi2)}){pair.second->resize(np*np);download(pair.first,pair.second->data());}
 o.chain=View("coordinate chain kernel cache",6*stride);o.positions=View("prolate position cache",2*stride);o.trig=View("azimuthal coordinate cache",2*np);
 auto ca=o.spectral.coordinate[0],cb=o.spectral.coordinate[1],chain=o.chain,position=o.positions,trig=o.trig;
 Kokkos::parallel_for("coordinate chain cache",Range(0,stride),KOKKOS_LAMBDA(int q){double c[6],x[2];prolate_coordinate_cache(ca(q%na),cb(q/na),b,c,x);for(int d=0;d<6;d++)chain(6*q+d)=c[d];position(2*q)=x[0];position(2*q+1)=x[1];});
 Kokkos::parallel_for("azimuthal coordinate cache",Range(0,np),KOKKOS_LAMBDA(int k){trig(2*k)=std::cos(2*Pi*k/np);trig(2*k+1)=std::sin(2*Pi*k/np);});
 o.spectral.mapping=View("regular factor kernel cache",15*points);
 auto ma=first[0],mma=second[0],mb=first[1],mmb=second[1],mapping=o.spectral.mapping;
 Kokkos::parallel_for("regular coefficient cache",Range(0,points),KOKKOS_LAMBDA(int p){
  int i=p%na,j=(p/na)%nb,mode=p/stride,m=mode<=np/2?mode:mode-np/2,r=hispid::AxisDerivatives::exponent(m);
  double a=.5*(ca(i)+1),B=cb(j),eta=-2*B/(1+B*B),den=1-eta*eta,S=hispid::AxisDerivatives::radial(a,r)*hispid::AxisDerivatives::angular(eta,r);
  double la=.5*(1/(1+a)+r/a),lla=-.25*(1/((1+a)*(1+a))+r/(a*a));
  double etaB=-2*(1-B*B)/std::pow(1+B*B,2),etaBB=4*B*(3-B*B)/std::pow(1+B*B,3),le=-r*eta/den,lle=-r*(1+eta*eta)/(den*den),lb=le*etaB,llb=lle*etaB*etaB+le*etaBB;
  int z=15*p;mapping(z)=S;mapping(z+1)=S*ma(i);mapping(z+2)=S*la;mapping(z+3)=S*mb(j);mapping(z+4)=S*lb;
  mapping(z+5)=S*ma(i)*ma(i);mapping(z+6)=S*(mma(i)+2*la*ma(i));mapping(z+7)=S*(la*la+lla);
  mapping(z+8)=S*mb(j)*mb(j);mapping(z+9)=S*(mmb(j)+2*lb*mb(j));mapping(z+10)=S*(lb*lb+llb);
  mapping(z+11)=S*ma(i)*mb(j);mapping(z+12)=S*la*mb(j);mapping(z+13)=S*lb*ma(i);mapping(z+14)=S*la*lb;
 });Exec().fence();
}
}
