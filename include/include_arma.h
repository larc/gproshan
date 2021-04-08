#ifndef INCLUDE_ARMA_H
#define INCLUDE_ARMA_H

#include "include.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

#ifdef NDEBUG
	#define ARMA_NO_DEBUG
#endif

// geometry processing and shape analysis framework
namespace gproshan {


#ifdef GPROSHAN_FLOAT
	typedef arma::fmat a_mat;
	typedef arma::fvec a_vec;
	typedef arma::frowvec a_rowvec;
	typedef arma::sp_fmat a_sp_mat;
#else
	typedef arma::mat a_mat;
	typedef arma::vec a_vec;
	typedef arma::rowvec a_rowvec;
	typedef arma::sp_mat a_sp_mat;
#endif


} // namespace gproshan

#endif // INCLUDE_ARMA_H

