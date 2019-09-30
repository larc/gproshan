#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "che.h"
#include "include_arma.h"

#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> sp_mat_e;

void laplacian(che * mesh, a_sp_mat & L, a_sp_mat & A);

void laplacian(che * mesh, sp_mat_e & L, sp_mat_e & A);

size_t eigs_laplacian(a_vec & eigval, a_mat & eigvec, che * mesh, const a_sp_mat & L, const a_sp_mat & A, const size_t & K);

#endif // LAPLACIAN_H

