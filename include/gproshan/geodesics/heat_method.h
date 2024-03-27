/*
 * Base on the code https://github.com/larc/dgpdec-course/tree/master/Geodesics
 * forked from https://github.com/dgpdec/course
 *
 * Geodesics in Heat: A New Approach to Computing Distance Based on Heat Method
 * Keenan Crane, Clarisse Weischedel, Max Wardetzky
 * To appear at ACM Transactions on Graphics
 */

#ifndef HEAT_METHOD_H
#define HEAT_METHOD_H

#include <gproshan/mesh/che.h>

#include <armadillo>


// geometry processing and shape analysis framework
namespace gproshan {

enum hm_opt {
				HEAT_ARMA,
				HEAT_CHOLMOD,
			#ifdef GPROSHAN_CUDA
				HEAT_CUDA
			#endif // GPROSHAN_CUDA
				};

double heat_method(float * dist, const che * mesh, const std::vector<index_t> & sources, const hm_opt & opt);

arma::vec compute_divergence(const che * mesh, const arma::vec & u);


#ifdef GPROSHAN_CUDA

double solve_positive_definite_gpu(arma::mat & x, const arma::sp_mat & A, const arma::mat & b);
double solve_positive_definite_cusolver(const int m, const int nnz, const double * hA_values, const int * hA_col_ptrs, const int * hA_row_indices, const double * hb, double * hx);

#endif // GPROSHAN_CUDA


} // namespace gproshan

#endif // HEAT_METHOD_H

