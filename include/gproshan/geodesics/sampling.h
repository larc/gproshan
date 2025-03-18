#ifndef SAMPLING_H
#define SAMPLING_H

#include <gproshan/geodesics/geodesics.h>
#include <gproshan/geodesics/geodesics_ptp.h>

#include <vector>


// geometry processing and shape analysis framework
namespace gproshan {


bool load_sampling(std::vector<index_t> & points, float & radio, che * mesh, size_t M);


} // namespace gproshan

#endif // SAMPLING_H

