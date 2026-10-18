#ifndef PUNCTURE_POLYNOMIAL_PRECONDITIONER_HPP
#define PUNCTURE_POLYNOMIAL_PRECONDITIONER_HPP
namespace puncture {
// Two stationary residual corrections of the same frozen linear system.
// The caller supplies backend vector operations and preallocated workspace.
// With two sweeps this is omega M^-1 (2I-omega A M^-1), a fixed linear
// right preconditioner. It changes neither A nor the nonlinear equations.
template<class Input,class Vector,class Inverse,class Action,class Copy,class Scale,class Axpy>
void relaxed_inverse(Input input,Vector output,Vector residual,Vector action,Vector correction,
                     Inverse inverse,Action apply,Copy copy,Scale scale,Axpy axpy,int sweeps,double omega){
 inverse(input,output);if(sweeps==1)return;scale(output,omega);
 for(int k=1;k<sweeps;k++){
  apply(output,action);copy(residual,input);axpy(residual,-1.,action);
  inverse(residual,correction);axpy(output,omega,correction);
 }
}
}
#endif
