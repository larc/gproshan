#include "d_basis_cosine.h"

#include <cassert>

// mesh dictionary learning and sparse coding namespace
namespace mdict {

basis_cosine::basis_cosine(const size_t & _r, const size_t & _n, const distance_t & _radio)
{
	radio = _radio;
	r = _r;
	n = _n;
	dim = r * n;
}

void basis_cosine::discrete(a_mat & phi, const a_mat & xy)
{
	assert(phi.n_cols == dim);

	a_vec x = xy.row(0).t();
	a_vec y = xy.row(1).t();

	vertex_t d = 1.0 / (r - 1);
	vertex_t c;

	for(size_t k = 0, ni = 1; ni <= n; ni++ )
	for(vertex_t alpha = 0; alpha <= 1; alpha += d, k++)
	{
		c = ni * M_PI / radio;
		phi.col(k) = cosine(x, y, c, alpha);
	}
}

void basis_cosine::plot_basis(ostream & os)
{
	vertex_t d = 1.0 / (r - 1);
	vertex_t c;

	os << "set multiplot layout " << n << "," << r << " rowsfirst scale 1.2;" << endl;

	for(size_t ni = 1; ni <= n; ni++ )
	for(vertex_t alpha = 0; alpha <= 1; alpha += d)
	{
		c = ni * M_PI / radio;
		os << "splot v * cos(u), v * sin(u), "; cosine(os, c, alpha); os << ";" << endl;
	}
}

void basis_cosine::plot_atoms(ostream & os, const a_vec & A)
{
	vertex_t d = 1.0 / (r - 1);
	vertex_t c;

	for(size_t k = 0, ni = 1; ni <= n; ni++ )
	for(vertex_t alpha = 0; alpha <= 1; alpha += d, k++)
	{
		c = ni * M_PI / radio;
		os << " + " << A(k) << " * "; cosine(os, c, alpha);
	}
}

a_vec basis_cosine::cosine(const a_vec & x, const a_vec & y, const vertex_t & c, const vertex_t & alpha)
{
	return cos(c * (alpha * x + (1 - alpha) * y));
}

void basis_cosine::cosine(ostream & os, const vertex_t & c, const vertex_t & alpha)
{
	os << "cos( " << c << " * (" << alpha << " * v * cos(u) + ( 1 - " << alpha << ") * v * sin(u)))";
}

} // mdict

