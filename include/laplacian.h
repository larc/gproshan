#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "che.h"

#include <armadillo>
#include <Eigen/Sparse>

using namespace arma;
using namespace Eigen;

typedef SparseMatrix<double> sp_mat_e;

void laplacian(che * mesh, sp_mat & L, sp_mat & A);

void laplacian(che * mesh, sp_mat_e & L, sp_mat_e & A);

void eigs_laplacian(vec & eigval, mat & eigvec, che * mesh, const sp_mat & L, const size_t & K);

#endif // LAPLACIAN_H

