/* Private-TU comparison against the untouched inherited GSL line functions. */
#include "../src/TP_Newton.c"
int main(void) {
  TwoPunctures_params_set_default();
  params_set_real("par_m_plus",.6); params_set_real("par_m_minus",.4);
  params_set_real("par_b",3); params_set_real("par_P_plus2",.03);
  params_set_real("par_S_plus3",.108); params_set_int("verbose",0);
  int checks=0;
  for(int fixture=0;fixture<2;fixture++) {
    int nv=1,na=fixture?10:8,nb=fixture?14:12,np=fixture?8:6,N=nv*na*nb*np;
    derivs *v,*u; allocate_derivs(&v,N); allocate_derivs(&u,N);
    double *F=dvector(0,N-1),*rhs=dvector(0,N-1),*a=dvector(0,N-1),*b=dvector(0,N-1);
    int *nc=ivector(0,N-1),**cols=imatrix(0,N-1,0,StencilSize-1);
    double **J=dmatrix(0,N-1,0,StencilSize-1);
    for(int refresh=0;refresh<2;refresh++) {
      for(int q=0;q<N;q++){v->d0[q]=(refresh+1)*1e-5*cos(.23*q);rhs[q]=sin(.37*q)+cos(.19*q);a[q]=b[q]=1e-4*sin(.13*q);}
      F_of_v(nv,na,nb,np,v,F,u); SetMatrix_JFD(nv,na,nb,np,u,nc,cols,J);
      TP_LineCache *cache=TP_cache_create(nv,na,nb,np,nc,cols,J); if(!cache)return 2;
      for(int repeat=0;repeat<2;repeat++) {
        for(int k=0;k<np;k++) {
          for(int i=0;i<na;i++) {
            LineRelax_be(a,i,k,nv,na,nb,np,rhs,nc,cols,J);
            TP_line_solve(cache->be+nv*(i+na*k),nv,b,rhs,nc,cols,J);
            if(memcmp(a,b,N*sizeof(double))){fprintf(stderr,"polar mismatch %d/%d/%d\n",fixture,refresh,i);return 3;}checks++;
          }
          for(int j=0;j<nb;j++) {
            LineRelax_al(a,j,k,nv,na,nb,np,rhs,nc,cols,J);
            TP_line_solve(cache->al+nv*(j+nb*k),nv,b,rhs,nc,cols,J);
            if(memcmp(a,b,N*sizeof(double))){fprintf(stderr,"radial mismatch\n");return 4;}checks++;
          }
        }
        for(int q=0;q<N;q++)rhs[q]=cos(.21*q+repeat+.7);
      }
      for(int it=0;it<NRELAX;it++) {
        relax(a,nv,na,nb,np,rhs,nc,cols,J,NULL); relax(b,nv,na,nb,np,rhs,nc,cols,J,cache);
        if(memcmp(a,b,N*sizeof(double))){fprintf(stderr,"sweep mismatch %d/%d/%d\n",fixture,refresh,it);return 5;}checks++;
      }
      TP_cache_destroy(cache);
      int saved=nc[0];nc[0]=65;cache=TP_cache_create(nv,na,nb,np,nc,cols,J);nc[0]=saved;
      if(cache){TP_cache_destroy(cache);return 7;}checks++;
      /* Force unsupported pivots: caller must receive NULL and use legacy GSL. */
      for(int m=0;m<nc[0];m++)if(cols[0][m]==0)J[0][m]=0;
      cache=TP_cache_create(nv,na,nb,np,nc,cols,J); if(cache){TP_cache_destroy(cache);return 6;}checks++;
    }
    free_derivs(v);free_derivs(u);free_dvector(F,0,N-1);free_dvector(rhs,0,N-1);free_dvector(a,0,N-1);free_dvector(b,0,N-1);
    free_ivector(nc,0,N-1);free_imatrix(cols,0,N-1,0,StencilSize-1);free_dmatrix(J,0,N-1,0,StencilSize-1);
  }
  printf("BY line/cache bitwise controls: %d passed (including 800 sweeps, refreshed JFD, unsupported-pivot/width cache rejection)\n",checks);
  return 0;
}
