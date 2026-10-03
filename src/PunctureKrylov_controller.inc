/* Shared C-compatible Krylov controller. Vector pointers belong to the
 * execution policy; this file never dereferences them. Small scalar banks
 * remain host allocations. Historical C arithmetic is in the reference policy. */
static int actual(int n,const double*b,const double*x,double*Ax,double*r,
                  PK_Apply A,void*c,PK_Result*out){
  if(!finite_vector(x,n))return PK_NONFINITE;
  int status=action(A,c,x,Ax,n,&out->operator_calls);if(status)return status;
  vector_difference(r,b,Ax,n);
  out->true_residual=length(r,n);double rhs=length(b,n);
  out->relative_residual=rhs>0?out->true_residual/rhs:out->true_residual;
  return isfinite(out->true_residual)?PK_SUCCESS:PK_NONFINITE;
}
static int gmres(int n,const double*b,double*x,const PK_Options*o,PK_Apply A,
                 PK_Apply M,void*c,PK_Monitor monitor,void*mc,PK_Result*out){
  int m=o->restart,status=PK_LIMIT;double *Ax=NULL,*r=NULL,*w=NULL,**V=NULL,**Z=NULL;
  double *H=NULL,*cs=NULL,*sn=NULL,*g=NULL,*y=NULL;
  vector_zero(x,n);double target=o->absolute_tolerance,beta=length(b,n);
  if(monitor)monitor(mc,0,beta,0,0,0);
  if(beta<=target){out->true_residual=beta;out->relative_residual=beta==0?0:1;out->recurrence_residual=beta;return PK_SUCCESS;}
  Ax=vector_bank(n);r=vector_bank(n);w=vector_bank(n);
  V=(double**)bank((size_t)m+1,1,sizeof(double*));Z=(double**)bank(m,1,sizeof(double*));H=(double*)bank((size_t)m+1,m,sizeof(double));
  cs=(double*)bank(m,1,sizeof(double));sn=(double*)bank(m,1,sizeof(double));g=(double*)bank((size_t)m+1,1,sizeof(double));y=(double*)bank(m,1,sizeof(double));
  if(!Ax||!r||!w||!V||!Z||!H||!cs||!sn||!g||!y){status=PK_ALLOCATION;goto done;}
  V[0]=vector_bank(n);if(!V[0]){status=PK_ALLOCATION;goto done;}
  if(o->eager_gmres_basis){
    for(int i=1;i<=m;i++)if(!(V[i]=vector_bank(n))){status=PK_ALLOCATION;goto done;}
    for(int i=0;i<m;i++)if(!(Z[i]=vector_bank(n))){status=PK_ALLOCATION;goto done;}
  }
  while(out->iterations<o->max_iterations){
    status=actual(n,b,x,Ax,r,A,c,out);if(status)goto done;beta=out->true_residual;
    if(beta<=target){status=PK_SUCCESS;goto done;}
    vector_divide(V[0],r,beta,n);
    memset(g,0,((size_t)m+1)*sizeof(double));g[0]=beta;memset(H,0,((size_t)m+1)*m*sizeof(double));
    int used=0;
    for(int k=0;k<m&&out->iterations<o->max_iterations;k++){
      if(!Z[k]&&!(Z[k]=vector_bank(n))){status=PK_ALLOCATION;goto done;}
      if(!V[k+1]&&!(V[k+1]=vector_bank(n))){status=PK_ALLOCATION;goto done;}
      status=action(M,c,V[k],Z[k],n,&out->preconditioner_calls);if(status)goto done;
      status=action(A,c,Z[k],w,n,&out->operator_calls);if(status)goto done;
      for(int pass=0;pass<2;pass++)for(int j=0;j<=k;j++){
        double a=dot(w,V[j],n);H[(size_t)j*m+k]+=a;vector_subtract(w,V[j],a,n);
      }
      H[(size_t)(k+1)*m+k]=length(w,n);
      if(H[(size_t)(k+1)*m+k]>0)vector_divide(V[k+1],w,H[(size_t)(k+1)*m+k],n);
      for(int j=0;j<k;j++){double a=cs[j]*H[(size_t)j*m+k]+sn[j]*H[(size_t)(j+1)*m+k];H[(size_t)(j+1)*m+k]=-sn[j]*H[(size_t)j*m+k]+cs[j]*H[(size_t)(j+1)*m+k];H[(size_t)j*m+k]=a;}
      double dd=hypot(H[(size_t)k*m+k],H[(size_t)(k+1)*m+k]);
      if(!(dd>0)||!isfinite(dd)){status=PK_BREAKDOWN;goto done;}
      cs[k]=H[(size_t)k*m+k]/dd;sn[k]=H[(size_t)(k+1)*m+k]/dd;H[(size_t)k*m+k]=dd;H[(size_t)(k+1)*m+k]=0;
      g[k+1]=-sn[k]*g[k];g[k]*=cs[k];used=k+1;out->iterations++;out->recurrence_residual=fabs(g[k+1]);
      if(monitor)monitor(mc,out->iterations,out->recurrence_residual,0,0,0);
      if(fabs(g[k+1])<=target)break;
    }
    for(int j=used-1;j>=0;j--){y[j]=g[j];for(int k=j+1;k<used;k++)y[j]-=H[(size_t)j*m+k]*y[k];y[j]/=H[(size_t)j*m+j];}
    for(int j=0;j<used;j++)vector_add(x,Z[j],y[j],n);
  }
  status=actual(n,b,x,Ax,r,A,c,out);if(!status)status=out->true_residual<=target?PK_SUCCESS:PK_LIMIT;
done:
  if(V)for(int i=0;i<=m;i++)vector_free(V[i]);if(Z)for(int i=0;i<m;i++)vector_free(Z[i]);
  vector_free(Ax);vector_free(r);vector_free(w);free(V);free(Z);free(H);free(cs);free(sn);free(g);free(y);return status;
}
static int bicgstab(int n,const double*b,double*x,const PK_Options*o,PK_Apply A,
                    PK_Apply M,void*c,PK_Monitor monitor,void*mc,PK_Result*out){
  double *work=vector_bank((size_t)8*n);if(!work)return PK_ALLOCATION;
  double *p=work,*rt=p+n,*s=rt+n,*t=s+n,*r=t+n,*v=r+n,*ph=v+n,*sh=ph+n;
  double alpha=0,beta=0,rho=0,rho1=1,omega=0;
  int status=action(A,c,x,r,n,&out->operator_calls);if(status)goto done;
  vector_shadow(rt,r,b,n);out->recurrence_residual=norm(o,r,n);
  if(monitor)monitor(mc,0,out->recurrence_residual,0,0,0);
  status=out->recurrence_residual<=o->absolute_tolerance?PK_SUCCESS:PK_LIMIT;
  if(out->recurrence_residual>o->absolute_tolerance)for(int ii=0;ii<o->max_iterations;ii++){
    out->iterations++;rho=inner(o,rt,r,n);
    if(o->legacy_bicgstab?fabs(rho)<1e-50:rho==0||!isfinite(rho)){status=PK_BREAKDOWN;break;}
    if(ii==0)vector_copy(p,r,n);
    else{beta=(rho/rho1)*(alpha/omega);vector_search(p,r,v,beta,omega,n);}
    vector_zero(ph,n);status=action(M,c,p,ph,n,&out->preconditioner_calls);if(status)break;
    status=action(A,c,ph,v,n,&out->operator_calls);if(status)break;
    double denominator=inner(o,rt,v,n);
    if(!o->legacy_bicgstab&&(denominator==0||!isfinite(denominator))){status=PK_BREAKDOWN;break;}
    alpha=rho/denominator;vector_step(s,r,v,alpha,n);out->recurrence_residual=norm(o,s,n);
    if(out->recurrence_residual<=o->absolute_tolerance){
      vector_add_left(x,ph,alpha,n);
      if(monitor)monitor(mc,ii+1,out->recurrence_residual,alpha,beta,omega);status=PK_SUCCESS;break;
    }
    vector_zero(sh,n);status=action(M,c,s,sh,n,&out->preconditioner_calls);if(status)break;
    status=action(A,c,sh,t,n,&out->operator_calls);if(status)break;
    denominator=inner(o,t,t,n);if(!o->legacy_bicgstab&&(denominator==0||!isfinite(denominator))){status=PK_BREAKDOWN;break;}
    omega=inner(o,t,s,n)/denominator;
    vector_bicg_update(x,r,s,t,ph,sh,alpha,omega,n);
    out->recurrence_residual=norm(o,r,n);
    if(monitor)monitor(mc,ii+1,out->recurrence_residual,alpha,beta,omega);
    if(out->recurrence_residual<=o->absolute_tolerance){status=PK_SUCCESS;break;}
    rho1=rho;
    if(o->legacy_bicgstab?fabs(omega)<1e-50:omega==0||!isfinite(omega)){status=PK_BREAKDOWN;break;}
    status=PK_LIMIT;
  }
  if(o->verify_true_residual){
    int verified=actual(n,b,x,v,r,A,c,out);
    if(!verified&&o->norm){
      out->true_residual=o->norm(r,n);double rhs=o->norm(b,n);
      out->relative_residual=rhs>0?out->true_residual/rhs:out->true_residual;
      if(!isfinite(out->true_residual))verified=PK_NONFINITE;
    }
    if(verified)status=verified;
    else if(status==PK_SUCCESS&&out->true_residual>o->absolute_tolerance)status=PK_LIMIT;
  }
done:vector_free(work);return status;
}
static int controller_solve(int n,const double*b,double*x,const PK_Options*o,PK_Apply A,
             PK_Apply M,void*c,PK_Monitor monitor,void*mc,PK_Result*out){
  if(!out)return PK_INVALID;memset(out,0,sizeof(*out));out->true_residual=out->relative_residual=NAN;
  if(n<=0||!b||!x||b==x||!o||!A||!M||o->method<0||o->method>1||o->max_iterations<1||o->restart<1||o->restart>4096||
     o->verify_true_residual<0||o->verify_true_residual>1||o->legacy_bicgstab<0||o->legacy_bicgstab>1||o->eager_gmres_basis<0||o->eager_gmres_basis>1||
     (o->method==PK_BICGSTAB&&!o->legacy_bicgstab&&!o->verify_true_residual)||
     !isfinite(o->absolute_tolerance)||o->absolute_tolerance<0||!finite_vector(b,n)||!finite_vector(x,n))return out->status=PK_INVALID;
  out->status=o->method==PK_GMRES?gmres(n,b,x,o,A,M,c,monitor,mc,out):bicgstab(n,b,x,o,A,M,c,monitor,mc,out);
  return out->status;
}
