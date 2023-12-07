#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include <gproshan/mesh/che.h>
#include <gproshan/include_arma.h>


// geometry processing and shape analysis framework
namespace gproshan {


void laplacian(const che * mesh, a_sp_mat & L, a_sp_mat & A);

size_t eigs_laplacian(const che * mesh, a_vec & eigval, a_mat & eigvec, a_sp_mat & L, a_sp_mat & A, const size_t k);


} // namespace gproshan

#endif // LAPLACIAN_H

