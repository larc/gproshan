#ifndef SAMPLING_H
#define SAMPLING_H

#include <gproshan/geodesics/geodesics.h>
#include <gproshan/geodesics/geodesics_ptp.h>

#include <vector>


// geometry processing and shape analysis framework
namespace gproshan {


index_t ** sampling_shape(std::vector<index_t> & points, size_t *& sizes, vertex *& normals, che * mesh, size_t n_points, real_t radio);

bool load_sampling(std::vector<index_t> & points, real_t & radio, che * mesh, size_t M);


} // namespace gproshan

#endif // SAMPLING_H

