#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "mesh/che.h"
#include "include_arma.h"

#include <Eigen/Sparse>


// geometry processing and shape analysis framework
namespace gproshan {


typedef Eigen::SparseMatrix<double> sp_mat_e;

void laplacian(che * mesh, a_sp_mat & L, a_sp_mat & A);

void laplacian(che * mesh, sp_mat_e & L, sp_mat_e & A);

size_t eigs_laplacian(a_vec & eigval, a_mat & eigvec, che * mesh, const a_sp_mat & L, const a_sp_mat & A, const size_t & K);


} // namespace gproshan

#endif // LAPLACIAN_H

