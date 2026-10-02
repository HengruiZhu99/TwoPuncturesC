/* Fixed-JFD line factors for the inherited TwoPunctures relaxation.
 * The Thomas elimination order follows GNU GSL 2.8 linalg/tridiag.c.
 * No reciprocal substitution, pivoting, or relaxation-order changes.
 * Source reference: https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz
 * GSL source is GPL-3.0-or-later; the recurrence here preserves its operation order.
 */
#ifndef TP_LINECACHE_H
#define TP_LINECACHE_H
#include <stdint.h>
typedef struct {
  int n, *index;
  uint64_t *off;
  double *alpha, *t, *above, *z, *x;
} TP_Line;
typedef struct { int nal, nbe; TP_Line *al, *be; } TP_LineCache;

static void TP_cache_destroy(TP_LineCache *c) {
  if (!c) return;
  for (int axis=0;axis<2;axis++) {
    TP_Line *bank=axis?c->be:c->al; int count=axis?c->nbe:c->nal;
    if (bank) for (int l=0;l<count;l++) {
      free(bank[l].index); free(bank[l].off); free(bank[l].alpha);
      free(bank[l].t); free(bank[l].above); free(bank[l].z); free(bank[l].x);
    }
    free(bank);
  }
  free(c);
}
static TP_LineCache *TP_cache_create(int nv,int na,int nb,int np,
                                    int *ncols,int **cols,double **J) {
  TP_LineCache *c=calloc(1,sizeof(*c)); if(!c)return NULL;
  c->nal=nv*nb*np; c->nbe=nv*na*np;
  c->al=calloc(c->nal,sizeof(TP_Line)); c->be=calloc(c->nbe,sizeof(TP_Line));
  if(!c->al||!c->be){TP_cache_destroy(c);return NULL;}
  for(int axis=0;axis<2;axis++) {
    int n=axis?nb:na, fixed=axis?na:nb;
    TP_Line *bank=axis?c->be:c->al;
    for(int k=0;k<np;k++)for(int f=0;f<fixed;f++)for(int v=0;v<nv;v++) {
      TP_Line *l=bank+v+nv*(f+fixed*k); l->n=n;
      l->index=malloc(n*sizeof(int)); l->off=calloc(n,sizeof(uint64_t));
      l->alpha=calloc(n,sizeof(double)); l->t=calloc(n,sizeof(double));
      l->above=calloc(n,sizeof(double)); l->z=malloc(n*sizeof(double)); l->x=malloc(n*sizeof(double));
      if(!l->index||!l->off||!l->alpha||!l->t||!l->above||!l->z||!l->x){TP_cache_destroy(c);return NULL;}
      for(int q=0;q<n;q++) {
        int Ic=Index(v,axis?f:q,axis?q:f,k,nv,na,nb,np);
        int Ip=Index(v,axis?f:q+1,axis?q+1:f,k,nv,na,nb,np);
        int Im=Index(v,axis?f:q-1,axis?q-1:f,k,nv,na,nb,np);
        double below=0; l->index[q]=Ic;
        if(ncols[Ic]>64){TP_cache_destroy(c);return NULL;}
        for(int m=0;m<ncols[Ic];m++) {
          int col=cols[Ic][m];
          if(col!=Ip&&col!=Ic&&col!=Im) l->off[q]|=UINT64_C(1)<<m;
          else {
            /* Independent assignments preserve reflected endpoints and overwrite semantics. */
            if(col==Im&&q>0) below=J[Ic][m];
            if(col==Ic) l->alpha[q]=J[Ic][m];
            if(col==Ip&&q<n-1) l->above[q]=J[Ic][m];
          }
        }
        if(q>0) {
          l->t[q]=below/l->alpha[q-1];
          l->alpha[q]=l->alpha[q]-l->t[q]*l->above[q-1];
        }
        /* Unsupported pivots retain the original GSL path and its error handling. */
        if(l->alpha[q]==0||!isfinite(l->alpha[q])){TP_cache_destroy(c);return NULL;}
      }
    }
  }
  return c;
}
static void TP_line_solve(TP_Line *line,int nv,double *dv,const double *rhs,
                          int *ncols,int **cols,double **J) {
  for(int v=0;v<nv;v++) {
    TP_Line *l=line+v; int n=l->n;
    for(int q=0;q<n;q++) {
      int Ic=l->index[q]; double b=rhs[Ic];
      for(int m=0;m<ncols[Ic];m++)if(l->off[q]&(UINT64_C(1)<<m)) b-=J[Ic][m]*dv[cols[Ic][m]];
      l->z[q]=b;
    }
    for(int q=1;q<n;q++)l->z[q]=l->z[q]-l->t[q]*l->z[q-1];
    l->x[n-1]=l->z[n-1]/l->alpha[n-1];
    for(int q=n-2;q>=0;q--)l->x[q]=(l->z[q]-l->above[q]*l->x[q+1])/l->alpha[q];
    for(int q=0;q<n;q++)dv[l->index[q]]=l->x[q];
  }
}
#endif
