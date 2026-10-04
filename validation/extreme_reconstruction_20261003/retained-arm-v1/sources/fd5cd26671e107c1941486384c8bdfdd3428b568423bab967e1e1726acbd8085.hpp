#ifndef HISPID_CARTESIAN_ORACLE_HPP
#define HISPID_CARTESIAN_ORACLE_HPP
// Independent continuous Cartesian expression for the complete retained
// polynomial, using puncture distances and all six Hessians. No legacy
// A/B transform or differentiation matrix is reused by this oracle.
inline Fields cartesian_polynomial(const HiSpID_Data&s,const double*x){
 using hispid::Jet;Jet X=Jet::variable(x[0],1),Y=Jet::variable(x[1],2),Z=Jet::variable(x[2],3);
 Jet rp=sqrt((X-Jet(s.b))*(X-Jet(s.b))+Y*Y+Z*Z),rm=sqrt((X+Jet(s.b))*(X+Jet(s.b))+Y*Y+Z*Z);
 Jet D=(rp+rm+Jet(2*s.b))/Jet(2),t=Jet(1)-Jet(2*s.b)/D,eta=(rm-rp)/Jet(2*s.b);
 const double lambda=hispid::AxisDerivatives::radial_stretch,kappa=hispid::AxisDerivatives::angular_stretch;
 Jet z=Jet(2)*t/(Jet(lambda)+Jet(1-lambda)*t)-Jet(1),angle=eta*Jet(std::tanh(kappa));
 Jet zeta=hispid::unary(angle,std::atanh(angle.v),1/(1-angle.v*angle.v),2*angle.v/std::pow(1-angle.v*angle.v,2))/Jet(kappa);
 const int na=s.local.n[0],nb=s.local.n[1],np=s.local.n[2],half=np/2;
 std::vector<Jet>Ta(na),Tb(nb);Ta[0]=Tb[0]=1;Ta[1]=z;Tb[1]=zeta;
 for(int i=2;i<na;i++)Ta[i]=Jet(2)*z*Ta[i-1]-Ta[i-2];
 for(int j=2;j<nb;j++)Tb[j]=Jet(2)*zeta*Tb[j-1]-Tb[j-2];
 Ta[0]=Tb[0]=.5;Jet result[4];
 for(int k=0;k<np;k++){
  int m=k<=half?k:k-half,r=hispid::AxisDerivatives::exponent(m);Jet re(1),im(0);
  for(int q=0;q<m;q++){Jet next=re*Y-im*Z;im=re*Z+im*Y;re=next;}
  Jet G=k<=half?re:im;if(r!=m)G=G*power(Y*Y+Z*Z,.5L*(r-m));
  Jet factor=Jet(-4*s.b*std::sqrt((m==0||m==half?1.:2.)/np))*G/power(D,r+1);
  for(int v=0;v<4;v++){
   Jet P;for(int j=0;j<nb;j++)for(int i=0;i<na;i++)P=P+Jet(s.coefficients[v+4*(i+na*(j+nb*k))])*Ta[i]*Tb[j];
   result[v]=result[v]+factor*P;
  }
 }
 Fields out{};for(int v=0;v<4;v++){
  out[v][0]=result[v].v;for(int d=1;d<4;d++)out[v][d]=result[v].d[d];
  for(int d=0;d<6;d++)out[v][4+d]=result[v].h[hi[d]][hj[d]];
 }return out;
}
#endif
