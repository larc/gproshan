#include "d_basis_cosine.h"

basis_cosine::basis_cosine(const vertex_t &_radio, const vertex_t & _r)
{
	radio = _radio;
	r = _r;
}

void basis_cosine::discrete(mat & phi, const mat & xy)
{
	vec x = xy.row(0).t();
	vec y = xy.row(1).t();
	
	size_t K = phi.n_cols;
	size_t n = K / r;
	vertex_t d = 1.0 / (r - 1);
	vertex_t c;

	for(size_t k = 0, ni = 1; ni <= n; ni++ )
	for(vertex_t alpha = 0; alpha <= 1; alpha += d, k++)
	{
		c = ni * M_PI / radio;
		phi.col(k) = cosine(x, y, c, alpha);	
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

