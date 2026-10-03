#include "TwoPunctures.h"
#include "TP_Modal.h"
#include <gsl/gsl_blas.h>
#include <stdint.h>
#include <limits.h>
struct TP_Modal {
  int na,nb,np,groups;
  double *lu,*transfer,*lower,*forward,*modal,*y,*column,*row_scale;
  size_t *permutation;
};
void TP_modal_destroy(TP_Modal *m){
  if(!m)return;
  free(m->lu);free(m->transfer);free(m->lower);free(m->forward);
  free(m->row_scale);free(m->modal);free(m->y);free(m->column);free(m->permutation);
  free(m);
}
int TP_modal_factors(const TP_Modal*m,const double**lu,const double**transfer,
                    const double**lower,const double**scale,const double**fourier,
                    const size_t**permutation){
 if(!m||!lu||!transfer||!lower||!scale||!fourier||!permutation)return -1;
 *lu=m->lu;*transfer=m->transfer;*lower=m->lower;*scale=m->row_scale;*fourier=m->forward;*permutation=m->permutation;return 0;
}
int TP_modal_shape(const TP_Modal*m,int*out){if(!m||!out)return -1;out[0]=m->na;out[1]=m->nb;out[2]=m->np;return 0;}
static size_t block_offset(const TP_Modal*m,int mode,int j){return ((size_t)mode*m->nb+j)*m->na*m->na;}
static int product(size_t a,size_t b,size_t*out){
  if(b&&a>SIZE_MAX/b)return -1;*out=a*b;return 0;
}
static TP_Modal *allocate_modal(int na,int nb,int np){
  if(na<2||nb<2||np<4||np%2)return NULL;
  size_t meridian,total,rows,blocks,fourier,bytes;
  if(product(na,nb,&meridian)||product(meridian,np,&total)||total>INT_MAX||
     product(meridian,np/2+1,&rows)||product(rows,na,&blocks)||
     product(np,np,&fourier)||fourier>INT_MAX||
     product(blocks,3*sizeof(double),&bytes))return NULL;
  /* Reject excessive/overflowing private allocations before touching input rows. */
  size_t extra[]={rows,rows,fourier,total,meridian,(size_t)na};
  for(size_t i=0;i<sizeof(extra)/sizeof(extra[0]);i++){
    size_t part;
    if(product(extra[i],i==1?sizeof(size_t):sizeof(double),&part)||bytes>SIZE_MAX-part)return NULL;
    bytes+=part;
  }
  if(bytes>((size_t)4<<30))return NULL;
  TP_Modal*m=calloc(1,sizeof(*m));if(!m)return NULL;
  m->na=na;m->nb=nb;m->np=np;m->groups=np/2+1;
  m->lu=calloc(blocks,sizeof(double));m->transfer=calloc(blocks,sizeof(double));m->lower=calloc(blocks,sizeof(double));
  m->row_scale=malloc(rows*sizeof(double));
  m->forward=malloc(fourier*sizeof(double));m->modal=malloc(total*sizeof(double));
  m->y=malloc(meridian*sizeof(double));m->column=malloc((size_t)na*sizeof(double));
  m->permutation=malloc(rows*sizeof(size_t));
  if(!m->lu||!m->transfer||!m->lower||!m->row_scale||!m->forward||!m->modal||!m->y||!m->column||!m->permutation){TP_modal_destroy(m);return NULL;}
  for(int mode=0;mode<np;mode++)for(int k=0;k<np;k++){
    int frequency=mode<=np/2?mode:mode-np/2;
    double normalization=(frequency==0||frequency==np/2)?1/sqrt((double)np):sqrt(2./np);
    m->forward[mode*np+k]=normalization*(mode<=np/2?cos(2*Pi*frequency*k/np):sin(2*Pi*frequency*k/np));
  }
  return m;
}
static TP_Modal *factor_modal(TP_Modal*m){
  int na=m->na,nb=m->nb;
  for(int mode=0;mode<m->groups;mode++)for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
    double scale=0;
    for(int q=0;q<na;q++){
      size_t at=block_offset(m,mode,j)+(size_t)i*na+q;
      if(!isfinite(m->lu[at])||!isfinite(m->lower[at])||!isfinite(m->transfer[at])){TP_modal_destroy(m);return NULL;}
      scale=fmax(scale,fmax(fabs(m->lu[at]),fmax(fabs(m->lower[at]),fabs(m->transfer[at]))));
    }
    if(!isfinite(scale)||scale==0){TP_modal_destroy(m);return NULL;}
    m->row_scale[i+na*(j+nb*mode)]=scale;
    for(int q=0;q<na;q++){
      size_t at=block_offset(m,mode,j)+(size_t)i*na+q;
      m->lu[at]/=scale;m->lower[at]/=scale;m->transfer[at]/=scale;
    }
  }
  for(int mode=0;mode<m->groups;mode++)for(int j=0;j<nb;j++){
    double *data=m->lu+block_offset(m,mode,j);
    gsl_matrix_view A=gsl_matrix_view_array(data,na,na);
    if(j){
      gsl_matrix_const_view L=gsl_matrix_const_view_array(m->lower+block_offset(m,mode,j),na,na);
      gsl_matrix_const_view T=gsl_matrix_const_view_array(m->transfer+block_offset(m,mode,j-1),na,na);
      if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1,&L.matrix,&T.matrix,1,&A.matrix)){TP_modal_destroy(m);return NULL;}
    }
    gsl_permutation p={(size_t)na,m->permutation+((size_t)mode*nb+j)*na};gsl_permutation_init(&p);int sign=0;
    for(size_t q=0;q<(size_t)na*na;q++)if(!isfinite(data[q])){TP_modal_destroy(m);return NULL;}
    gsl_error_handler_t *handler=gsl_set_error_handler_off();
    int status=gsl_linalg_LU_decomp(&A.matrix,&p,&sign);gsl_set_error_handler(handler);
    if(status){TP_modal_destroy(m);return NULL;}
    for(int i=0;i<na;i++)if(!isfinite(data[i*na+i])||data[i*na+i]==0){TP_modal_destroy(m);return NULL;}
    if(j+1<nb)for(int q=0;q<na;q++){
      for(int i=0;i<na;i++)m->column[i]=m->transfer[block_offset(m,mode,j)+(size_t)i*na+q];
      gsl_vector_view b=gsl_vector_view_array(m->column,na);
      gsl_vector_view x=gsl_vector_view_array_with_stride(m->transfer+block_offset(m,mode,j)+q,na,na);
      handler=gsl_set_error_handler_off();status=gsl_linalg_LU_solve(&A.matrix,&p,&b.vector,&x.vector);gsl_set_error_handler(handler);
      if(status){TP_modal_destroy(m);return NULL;}
      for(int i=0;i<na;i++)if(!isfinite(m->transfer[block_offset(m,mode,j)+(size_t)i*na+q])){TP_modal_destroy(m);return NULL;}
    }
  }
  return m;
}
TP_Modal *TP_modal_create(int nv,int na,int nb,int np,
                          const int *ncols,int *const *cols,double *const *J){
  if(nv!=1||!ncols||!cols||!J)return NULL;
  TP_Modal*m=allocate_modal(na,nb,np);if(!m)return NULL;
  size_t total=(size_t)na*nb*np;
  for(size_t row=0;row<total;row++)if(ncols[row]<0||ncols[row]>StencilSize||!cols[row]||!J[row]){TP_modal_destroy(m);return NULL;}
  for(int mode=0;mode<m->groups;mode++)for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
    for(int k=0;k<np;k++){
      int row=i+na*(j+nb*k);
      if(ncols[row]<0||ncols[row]>StencilSize){TP_modal_destroy(m);return NULL;}
      for(int q=0;q<ncols[row];q++){
        int col=cols[row][q];
        if(col<0||(size_t)col>=total||!isfinite(J[row][q])){TP_modal_destroy(m);return NULL;}
        int ci=col%na,cj=(col/na)%nb,ck=col/(na*nb);
        double angle=2*Pi*mode*(ck-k)/np,value=J[row][q]/np;
        size_t at=block_offset(m,mode,j)+(size_t)i*na+ci;
        if(cj==j)m->lu[at]+=value*cos(angle);
        else if(cj==j-1)m->lower[at]+=value*cos(angle);
        else if(cj==j+1)m->transfer[at]+=value*cos(angle);
        else{TP_modal_destroy(m);return NULL;}
      }
    }
    /* cos projects onto the EVEN cyclic azimuthal average. This explicitly
     * symmetrizes opposite phi shifts; full J remains unchanged. */
  }
  return factor_modal(m);
}
/* Orthogonal prolate Laplacian applied to u=(A-1)V, discretized in
 * alpha/beta/phi exactly as the legacy JFD. Cross derivatives cancel
 * analytically. Boundary folding follows TP Index (no Fourier parity).
 * The full spectral residual/JVP still uses the original implementation. */
TP_Modal *TP_modal_create_analytic(int na,int nb,int np,const double*U){
  if(!U)return NULL;
  TP_Modal*m=allocate_modal(na,nb,np);if(!m)return NULL;
  double b=params_get_real("par_b"),mp=params_get_real("par_m_plus"),mm=params_get_real("par_m_minus");
  if(!isfinite(b)||b<=0||!isfinite(mp)||!isfinite(mm)){TP_modal_destroy(m);return NULL;}
  double ha=Pi/na,hb=Pi/nb,hp=2*Pi/np;
  for(int j=0;j<nb;j++)for(int i=0;i<na;i++){
    double al=ha*(i+.5),be=hb*(j+.5),sa=sin(al),sb=sin(be);
    double A=-cos(al),B=-cos(be),a=.5*(A+1),d=A-1;
    double X=2*atanh(a),R=Pih+2*atan(B),sx=sinh(X),cx=cosh(X),sr=sin(R),cr=cos(R);
    double H=b*b*(sx*sx+sr*sr),rho=b*sx*sr,x=b*cx*cr;
    double h=1-a*a,hX=-a*h,k=.5*(1+B*B),kR=B*k;
    double cAA=d*h*h/H,cA=(2*h*h+d*(hX+h*cx/sx))/H;
    double cBB=d*k*k/H,cB=d*(kR+k*cr/sr)/H;
    double aa=cAA/(sa*sa),ab=cA/sa-cAA*cos(al)/(sa*sa*sa);
    double bb=cBB/(sb*sb),bc=cB/sb-cBB*cos(be)/(sb*sb*sb);
    double w=pow(sa*sb,3),am=w*(aa/(ha*ha)-ab/(2*ha)),ap=w*(aa/(ha*ha)+ab/(2*ha));
    double bm=w*(bb/(hb*hb)-bc/(2*hb)),bp=w*(bb/(hb*hb)+bc/(2*hb));
    double potential=0;
    for(int phi=0;phi<np;phi++){
      double y=rho*cos(hp*phi),z=rho*sin(hp*phi);
      double rp=sqrt((x-b)*(x-b)+y*y+z*z),rm=sqrt((x+b)*(x+b)+y*y+z*z);
      double psi=1+.5*mp/rp+.5*mm/rm+U[i+na*(j+nb*phi)];
      if(!isfinite(psi)||psi<=0){TP_modal_destroy(m);return NULL;}
      double p2=psi*psi,p4=p2*p2;
      potential+=.875*BY_KKofxyz(x,y,z)/(p4*p4)/np;
    }
    double c0=(hX+h*cx/sx)/H-d*potential;
    for(int mode=0;mode<m->groups;mode++){
      double lambda=-4*pow(sin(Pi*mode/np),2)/(hp*hp);
      size_t at=block_offset(m,mode,j)+(size_t)i*na;
      m->lu[at+i]=w*c0-am-ap-bm-bp+w*d/(rho*rho)*lambda;
      if(i)m->lu[at+i-1]=am;else m->lu[at+i]+=am;
      if(i+1<na)m->lu[at+i+1]=ap;else m->lu[at+i]+=ap;
      if(j)m->lower[at+i]=bm;else m->lu[at+i]+=bm;
      if(j+1<nb)m->transfer[at+i]=bp;else m->lu[at+i]+=bp;
    }
  }
  return factor_modal(m);
}
static void forward(const TP_Modal*m,const double*b,double*out){
  int stride=m->na*m->nb,N=m->np;
  for(int row=0;row<stride;row++){
    double sum=0,error=0;
    for(int k=0;k<N;k++){double v=b[row+stride*k]-error,t=sum+v;error=(t-sum)-v;sum=t;}
    double mean=sum/N;
    for(int mode=0;mode<N;mode++){
      double value=mode==0?sum/sqrt((double)N):0;
      if(mode)for(int k=0;k<N;k++)value+=m->forward[mode*N+k]*(b[row+stride*k]-mean);
      out[row+stride*mode]=value;
    }
  }
}
static void backward(const TP_Modal*m,const double*b,double*out){
  int stride=m->na*m->nb,N=m->np;
  for(int row=0;row<stride;row++)for(int k=0;k<N;k++){
    double value=0;for(int mode=0;mode<N;mode++)value+=m->forward[mode*N+k]*b[row+stride*mode];out[row+stride*k]=value;
  }
}
int TP_modal_solve(TP_Modal*m,const double*b,double*x){
  if(!m||!b||!x)return -1;
  for(int q=0;q<m->na*m->nb*m->np;q++)if(!isfinite(b[q]))return -1;
  int na=m->na,nb=m->nb,stride=na*nb;forward(m,b,m->modal);
  for(int mode=0;mode<m->np;mode++){
    int group=mode<=m->np/2?mode:mode-m->np/2;
    for(int j=0;j<nb;j++){
      for(int i=0;i<na;i++){
        double value=m->modal[i+na*j+stride*mode]/m->row_scale[i+na*(j+nb*group)];
        if(j)for(int q=0;q<na;q++)value-=m->lower[block_offset(m,group,j)+(size_t)i*na+q]*m->y[q+na*(j-1)];
        m->column[i]=value;
      }
      gsl_matrix_const_view A=gsl_matrix_const_view_array(m->lu+block_offset(m,group,j),na,na);
      gsl_permutation p={(size_t)na,m->permutation+((size_t)group*nb+j)*na};
      gsl_vector_view rhs=gsl_vector_view_array(m->column,na),out=gsl_vector_view_array(m->y+na*j,na);
      gsl_error_handler_t *handler=gsl_set_error_handler_off();
      int status=gsl_linalg_LU_solve(&A.matrix,&p,&rhs.vector,&out.vector);gsl_set_error_handler(handler);
      if(status)return -1;
    }
    for(int j=nb-1;j>=0;j--)for(int i=0;i<na;i++){
      double value=m->y[i+na*j];
      if(j+1<nb)for(int q=0;q<na;q++)value-=m->transfer[block_offset(m,group,j)+(size_t)i*na+q]*m->modal[q+na*(j+1)+stride*mode];
      m->modal[i+na*j+stride*mode]=value;
    }
  }
  backward(m,m->modal,x);
  for(int q=0;q<m->na*m->nb*m->np;q++)if(!isfinite(x[q]))return -1;
  return 0;
}
