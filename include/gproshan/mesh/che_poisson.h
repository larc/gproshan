#ifndef CHE_POISSON_H
#define CHE_POISSON_H

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


/**
	Solve (-1)^k L(AL)^(k-1) X = (-1)^k A^(-1) B
*/
void poisson(che * mesh, const size_t old_n_vertices, index_t k);

void biharmonic_interp_2(che * mesh, const size_t old_n_vertices, const size_t n_vertices, const std::vector<index_t> & border_vertices, const index_t k);


} // namespace gproshan

#endif //CHE_POISSON_H

