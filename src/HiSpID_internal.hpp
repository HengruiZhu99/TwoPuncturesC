#ifndef HISPID_INTERNAL_HPP
#define HISPID_INTERNAL_HPP
#include "HiSpID.h"
#include "HiSpID_jets.hpp"
#include <array>
#include <vector>
#include <string>
namespace hispid {
struct Seed { Jet metric[3][3], physical[3][3], extrinsic[3][3], A[3][3], psi,K; };
struct Background {
 Jet metric[3][3],inv[3][3],C[3][3][3],M[3][3],psi,K,g,far_correction;
 Jet opmetric[3][3],opinv[3][3],opC[3][3][3];
 double R,lapPsi,divM[3];
};
void seed(const HiSpID_Hole&,int,const double*,Seed&);
void background(const HiSpID_Config&,const double*,Background&);
bool valid(const HiSpID_Config&);
void longitudinal(const Jet metric[3][3],const Jet C[3][3][3],
                   const Jet b[3],Jet L[3][3]);
void divergence(const Jet inv[3][3],const Jet C[3][3][3],
                const Jet T[3][3],double out[3]);
double laplacian(const Jet inv[3][3],const Jet C[3][3][3],const Jet&u);
void equations(const Background&,const Jet fields[4],double out[4],
               const Jet *direction=nullptr);
extern thread_local std::string last_error;
}
#endif
