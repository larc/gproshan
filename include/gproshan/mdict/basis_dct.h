#ifndef BASIS_DCT_H
#define BASIS_DCT_H

#include <gproshan/mdict/basis.h>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class basis_dct: public basis
{
	private:
		size_t n_freq;		///< frequency

	public:
		basis_dct(const size_t n, const float r = 1);
		void discrete(arma::fmat & phi, const arma::fvec & x, const arma::fvec & y);
		void d_discrete(arma::fmat & phi, const arma::fvec & x, const arma::fvec & y, const bool b);
		float freq(const index_t idx);

	private:
		void plot_basis(std::ostream & os);
		void plot_atoms(std::ostream & os, const arma::fvec & A);
		arma::fvec dct(const arma::fvec & x, const arma::fvec & y, const index_t nx, const index_t ny);
		arma::fvec d_dct(const arma::fvec & x, const arma::fvec & y, const index_t nx, const index_t ny);
		void dct(std::ostream & os, const index_t nx, const index_t ny);
};


} // namespace gproshan::mdict

#endif // BASIS_DCT_H

