#ifndef BASIS_COSSINE_H
#define BASIS_COSSINE_H

#include <gproshan/mdict/basis.h>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class basis_cosine: public basis
{
	private:
		size_t n_rot;		///< rotations
		size_t n_freq;		///< frequency

	public:
		basis_cosine(const size_t nr, const size_t nf, const float r = 0);
		void discrete(arma::fmat & phi, const arma::fvec & x, const arma::fvec & y);

	private:
		void plot_basis(std::ostream & os);
		void plot_atoms(std::ostream & os, const arma::fvec & A);
		arma::fvec cosine(const arma::fvec & x, const arma::fvec & y, const float c, const float alpha);
		void cosine(std::ostream & os, const float c, const float alpha);
};


} // namespace gproshan::mdict

#endif // BASIS_COSSINE_H

