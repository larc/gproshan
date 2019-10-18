/*
 * Base on the code https://github.com/larc/dgpdec-course/tree/master/Geodesics
 * forked from https://github.com/dgpdec/course
 *
 * Geodesics in Heat: A New Approach to Computing Distance Based on Heat Flow
 * Keenan Crane, Clarisse Weischedel, Max Wardetzky
 * To appear at ACM Transactions on Graphics
 */

#ifndef HEAT_FLOW_H
#define HEAT_FLOW_H

#include "che.h"
#include "include_arma.h"

#include <cholmod.h>	// suitesparse/cholmod.h


// geometry processing and shape analysis framework
namespace gproshan {


distance_t * heat_flow(che * mesh, const std::vector<index_t> & sources, double & solve_time);

distance_t * heat_flow_gpu(che * mesh, const std::vector<index_t> & sources, double & solve_time);

void compute_divergence(che * mesh, const a_mat & u, a_mat & div);

/// cholmod Keenan implementation
/// base on the code https://github.com/larc/dgpdec-course/tree/master/Geodesics
double solve_positive_definite(a_mat & x, const a_sp_mat & A, const a_mat & b, cholmod_common * context);

cholmod_dense * arma_2_cholmod(const a_mat & m, cholmod_common * context);

cholmod_sparse * arma_2_cholmod(const a_sp_mat & m, cholmod_common * context);

/// 
double solve_positive_definite_gpu(a_mat & x, const a_sp_mat & A, const a_mat & b);

#ifdef CUDA_SUPPORT

/// host and device support
/// https://docs.nvidia.com/cuda/cusolver/index.html#cusolver-lt-t-gt-csrlsvchol
double solve_positive_definite_cusolver(const int m, const int nnz, const real_t * hA_values, const int * hA_col_ptrs, const int * hA_row_indices, const real_t * hb, real_t * hx, const bool host = 0);

/// device only, incomplete cholesky factorization
/// https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csric02
double solve_positive_definite_cusparse(const int m, const int nnz, const real_t * hA_values, const int * hA_col_ptrs, const int * hA_row_indices, const real_t * hb, real_t * hx);

/// host and device support
/// using cusolverSp_LOWLEVEL_PREVIEW.h library
/// no documentation, code base on cuda/samples/7_CUDALibraries/cuSolverSp_LowlevelCholesky
double solve_positive_definite_cusolver_preview(const int m, const int nnz, const real_t * hA_values, const int * hA_col_ptrs, const int * hA_row_indices, const real_t * hb, real_t * hx, const bool host = 0);

#endif


} // namespace gproshan

#endif // HEAT_FLOW_H

