#ifndef HISPID_INTERNAL_HPP
#define HISPID_INTERNAL_HPP
#include "HiSpID.h"
#include "HiSpID_geometry_kernels.hpp"
#include <array>
#include <vector>
#include <string>
namespace hispid {
void seed(const HiSpID_Hole&,int,const double*,Seed&);
void background(const HiSpID_Config&,const double*,Background&,int family=HISPID_SEED_QI);
double allocation_bound(const HiSpID_Config&,bool sampler_only=false);
bool valid(const HiSpID_Config&,bool sampler_only=false);
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
