#ifndef PUNCTURE_LEFT_GMRES_HPP
#define PUNCTURE_LEFT_GMRES_HPP
#include "PunctureKrylov.h"
#include <algorithm>
#include <cmath>
namespace puncture {
// Restarted left-preconditioned inner cycles. Success always requires an
// explicit residual in the original system, with the original norm/target.
// All workspace and vector primitives are supplied by the execution policy.
template<class V,class Length,class Difference,class Add,class Action,class Inverse,class Cycle,class Progress>
int left_gmres(V rhs,V x,V residual,V conditioned,V correction,V action,
               PK_Options options,double rtol,Length length,Difference difference,Add add,
               Action apply,Inverse inverse,Cycle cycle,Progress progress,PK_Result&result){
 const double original=length(rhs),target=rtol*original;
 result=PK_Result{};result.status=PK_LIMIT;
 apply(x,action);difference(residual,rhs,action);
 while(result.iterations<options.max_iterations){
  const double norm=length(residual);result.true_residual=norm;result.relative_residual=original?norm/original:norm;
  progress(result.iterations,result.relative_residual);
  if(!std::isfinite(norm)){result.status=PK_NONFINITE;break;}
  if(norm<=target){result.status=PK_SUCCESS;break;}
  inverse(residual,conditioned);const double scale=length(conditioned);
  if(!std::isfinite(scale)||scale==0){result.status=PK_BREAKDOWN;break;}
  PK_Options inner=options;inner.method=PK_GMRES;
  inner.max_iterations=std::min(options.restart,options.max_iterations-result.iterations);
  inner.absolute_tolerance=rtol*scale;
  PK_Result part{};
  int status=cycle(conditioned,correction,inner,[&](V in,V out){apply(in,action);inverse(action,out);},part,scale);
  result.iterations+=part.iterations;result.operator_calls+=part.operator_calls;
  result.preconditioner_calls+=part.preconditioner_calls;result.recurrence_residual=part.recurrence_residual;
  add(x,correction);apply(x,action);difference(residual,rhs,action);
  result.true_residual=length(residual);result.relative_residual=original?result.true_residual/original:result.true_residual;
  if(std::isfinite(result.true_residual)&&result.true_residual<=target){result.status=PK_SUCCESS;break;}
  if(status!=PK_SUCCESS&&status!=PK_LIMIT){result.status=status;break;}
  if(!part.iterations){result.status=PK_BREAKDOWN;break;}
 }
 progress(result.iterations,result.relative_residual);return result.status;
}
}
#endif
