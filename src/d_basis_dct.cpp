#include "d_basis_dct.h"

#include <cassert>

basis_dct::basis_dct(const vertex_t &_radio, const size_t & _n)
{
	radio = _radio;
	n = _n;
}

void basis_dct::discrete(mat & phi, const mat & xy)
{
	assert(phi.n_cols != n * n);

	vec x = xy.row(0).t();
	vec y = xy.row(1).t();
	
	for(index_t k = 0, nx = 0; nx < n; nx++)
	for(index_t ny = 0; ny < n; ny++, k++)
		phi.col(k) = dct(x, y, nx, ny);	
}

void basis_dct::plot_basis(ostream & os)
{
	os << "set multiplot layout " << n  << "," << n << " rowsfirst scale 1.2;" << endl;
	
	for(index_t nx = 0; nx < n; nx++)
	for(index_t ny = 0; ny < n; ny++)
	{
		os << "splot v * cos(u), v * sin(u), "; dct(os, nx, ny); os << ";" << endl;
	}
}

void basis_dct::plot_atoms(ostream & os, const vec & A)
{
	for(index_t k = 0, nx = 0; nx < n; nx++)
	for(index_t ny = 0; ny < n; ny++, k++)
	{
		os << " + " << A(k) << " * "; dct(os, nx, ny);
	}
}

vec basis_dct::dct(const vec & x, const vec & y, const index_t & nx, const index_t & ny)
{
	return cos( (M_PI * x * nx) / radio ) % cos( (M_PI * y * ny) / radio );
}

void basis_dct::dct(ostream & os, const index_t & nx, const index_t & ny)
{
	os << "cos( (pi * v * cos(u) * " << nx << " ) / " << radio << " ) *";
	os << "cos( (pi * v * sin(u) * " << ny << " ) / " << radio << " )";
}

