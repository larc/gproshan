#include "d_basis_cosine.h"

#include <cassert>

basis_cosine::basis_cosine(const vertex_t &_radio, const size_t & _r, const size_t & _n)
{
	radio = _radio;
	r = _r;
	n = _n;
	dim = r * n;
}

void basis_cosine::discrete(mat & phi, const mat & xy)
{
	assert(phi.n_cols == dim);

	vec x = xy.row(0).t();
	vec y = xy.row(1).t();
	
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
	
	os << "set multiplot layout " << n  << "," << r << " rowsfirst scale 1.2;" << endl;
	
	for(size_t ni = 1; ni <= n; ni++ )
	for(vertex_t alpha = 0; alpha <= 1; alpha += d)
	{
		c = ni * M_PI / radio;
		os << "splot v * cos(u), v * sin(u), "; cosine(os, c, alpha); os << ";" << endl;
	}
}

void basis_cosine::plot_atoms(ostream & os, const vec & A)
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

vec basis_cosine::cosine(const vec & x, const vec & y, const vertex_t & c, const vertex_t & alpha)
{
	return cos(c * (alpha * x + (1 - alpha) * y));
}

void basis_cosine::cosine(ostream & os, const vertex_t & c, const vertex_t & alpha)
{
	os << "cos( " << c << " * (" << alpha << " * v * cos(u) + ( 1 - " << alpha << ") * v * sin(u)))";
}

