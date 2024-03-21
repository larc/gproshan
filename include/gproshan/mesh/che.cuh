#ifndef CHE_CUH
#define CHE_CUH

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


void cuda_create_CHE(const che * h_che, che *& dd_che, che *& d_che, const bool normal = false, const bool color = false);

void cuda_free_CHE(che *& dd_che, che *& d_che);


} // namespace gproshan

#endif // CHE_CUH

