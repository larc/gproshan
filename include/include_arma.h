#include "include.h"

#include <armadillo>

#ifdef SINGLE_P
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

