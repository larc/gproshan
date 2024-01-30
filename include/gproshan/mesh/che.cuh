#ifndef CHE_CUH
#define CHE_CUH

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


void cuda_create_CHE(CHE * h_che, CHE *& dd_che, CHE *& d_che, const bool normal = false, const bool color = false);

void cuda_free_CHE(CHE *& dd_che, CHE *& d_che);


} // namespace gproshan

#endif // CHE_CUH

