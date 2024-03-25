#include <gproshan/mdict/basis_cosine.h>

#include <cassert>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


basis_cosine::basis_cosine(const size_t nr, const size_t nf, const real_t r): basis(r, r * nf), n_rot(nr), n_freq(nf) {}

void basis_cosine::discrete(arma::fmat & phi, const arma::fvec & x, const arma::fvec & y)
{
	assert(phi.n_cols == _dim);

	real_t d = 1.0 / (n_rot - 1);
	real_t c;

	for(size_t k = 0, ni = 1; ni <= n_freq; ++ni)
	for(real_t alpha = 0; alpha <= 1; alpha += d, ++k)
	{
		c = ni * M_PI / _radio;
		phi.col(k) = cosine(x, y, c, alpha);
	}
}

void basis_cosine::plot_basis(std::ostream & os)
{
	real_t d = 1.0 / (n_rot - 1);
	real_t c;

	os << "set multiplot layout " << n_freq << "," << n_rot << " rowsfirst scale 1.2;" << std::endl;

	for(size_t ni = 1; ni <= n_freq; ++ni)
	for(real_t alpha = 0; alpha <= 1; alpha += d)
	{
		c = ni * M_PI / _radio;
		os << "splot v * cos(u), v * sin(u), "; cosine(os, c, alpha); os << ";" << std::endl;
	}
}

void basis_cosine::plot_atoms(std::ostream & os, const arma::fvec & A)
{
	real_t d = 1.0 / (n_rot - 1);
	real_t c;

	for(size_t k = 0, ni = 1; ni <= n_freq; ++ni)
	for(real_t alpha = 0; alpha <= 1; alpha += d, ++k)
	{
		c = ni * M_PI / _radio;
		os << " + " << A(k) << " * "; cosine(os, c, alpha);
	}
}

arma::fvec basis_cosine::cosine(const arma::fvec & x, const arma::fvec & y, const real_t c, const real_t alpha)
{
	return cos(c * (alpha * x + (1 - alpha) * y));
}

void basis_cosine::cosine(std::ostream & os, const real_t c, const real_t alpha)
{
	os << "cos( " << c << " * (" << alpha << " * v * cos(u) + ( 1 - " << alpha << ") * v * sin(u)))";
}


} // namespace gproshan::mdict

