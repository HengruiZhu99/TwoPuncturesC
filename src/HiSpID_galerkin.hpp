#ifndef HISPID_GALERKIN_HPP
#define HISPID_GALERKIN_HPP
#include "PunctureModalProjection.hpp"
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>
namespace hispid {
// Cardinal interpolation at Chebyshev roots and its quadrature-weighted
// transpose. R=(n/q) B^T is the discrete Chebyshev L2 projection; R B=I.
struct GalerkinTransfer {
 std::array<int,3> n,q;
 int nv;
 std::array<std::vector<double>,2>B,R;
 mutable std::vector<double>scratch;
 GalerkinTransfer(const int*coarse,const int*fine,int variables=4):nv(variables){
  std::copy_n(coarse,3,n.begin());std::copy_n(fine,3,q.begin());
  if(n[2]!=q[2]||n[0]>q[0]||n[1]>q[1]||nv<1)throw std::invalid_argument("invalid Galerkin transfer shape");
  const long double pi=std::acos(-1.L);
  for(int axis=0;axis<2;axis++){
   int a=n[axis],b=q[axis];B[axis].resize(size_t(a)*b);R[axis].resize(size_t(a)*b);
   std::vector<long double>node(a),w(a),terms(a);
   for(int i=0;i<a;i++){long double t=pi*(i+.5L)/a;node[i]=-std::cos(t);w[i]=(i%2?-1.L:1.L)*std::sin(t);}
   for(int j=0;j<b;j++){
    long double x=-std::cos(pi*(j+.5L)/b),sum=0;int hit=-1;
    for(int i=0;i<a;i++)if(std::abs(x-node[i])<8*std::numeric_limits<long double>::epsilon()){hit=i;break;}
    if(hit>=0){for(int i=0;i<a;i++)B[axis][j*a+i]=i==hit?1:0;}
    else{
     for(int i=0;i<a;i++){terms[i]=w[i]/(x-node[i]);sum+=terms[i];}
     for(int i=0;i<a;i++)B[axis][j*a+i]=double(terms[i]/sum);
    }
    for(int i=0;i<a;i++)R[axis][i*b+j]=double(a)/b*B[axis][j*a+i];
   }
  }
 }
 static void along(int axis,int rows,int cols,int radial,int polar,int phi,int variables,
                   const std::vector<double>&matrix,const double*in,double*out){
  int ia=axis==0?cols:radial,ib=axis==0?polar:cols;
  int oa=axis==0?rows:radial,ob=axis==0?polar:rows;
  for(int k=0;k<phi;k++)for(int j=0;j<ob;j++)for(int i=0;i<oa;i++){
   int row=axis==0?i:j,anchor=std::min(cols-1,int((row+.5)*cols/rows));
   int base=axis==0?ia*(j+ib*k):i+ia*ib*k,stride=axis==0?1:ia;
   for(int v=0;v<variables;v++){
    const double*x=in+variables*base+v;double center=x[variables*stride*anchor];
    out[variables*(i+oa*(j+ob*k))+v]=center+puncture::compensated_projection(x,variables*stride,matrix.data()+row*cols,cols,center);
   }
  }
 }
 void prolong(const double*in,double*out)const{
  scratch.resize(size_t(nv)*q[0]*q[1]*q[2]);
  along(0,q[0],n[0],n[0],n[1],n[2],nv,B[0],in,scratch.data());
  along(1,q[1],n[1],q[0],n[1],n[2],nv,B[1],scratch.data(),out);
 }
 void restrict(const double*in,double*out)const{
  scratch.resize(size_t(nv)*q[0]*q[1]*q[2]);
  along(1,n[1],q[1],q[0],q[1],q[2],nv,R[1],in,scratch.data());
  along(0,n[0],q[0],q[0],n[1],q[2],nv,R[0],scratch.data(),out);
 }
 long double matrix_bytes()const{return 16.L*(n[0]*q[0]+n[1]*q[1]);}
};
}
#endif
