#ifndef SAMPLING_H
#define SAMPLING_H

#include "geodesics.h"
#include "geodesics_ptp.h"

#include <vector>


// geometry processing and shape analysis framework
namespace gproshan {


index_t ** sampling_shape(std::vector<index_t> & points, size_t *& sizes, vertex *& normals, che * shape, size_t n_points, distance_t radio);

bool load_sampling(std::vector<index_t> & points, distance_t & radio, che * mesh, size_t M);


} // namespace gproshan

#endif // SAMPLING_H

