#include "basis_dct.h"

#include <cassert>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


basis_dct::basis_dct(const size_t & _n, const distance_t & _radio)
{
	radio = _radio;
	n = _n;
	dim = n * n;
}

void basis_dct::discrete(a_mat & phi, const a_mat & xy)
{
	assert(phi.n_cols == dim);

	a_vec x = xy.row(0).t();
	a_vec y = xy.row(1).t();
	for(index_t k = 0, nx = 0; nx < n; nx++)
	for(index_t ny = 0; ny < n; ny++, k++)
	{
		phi.col(k) = dct(x, y, nx, ny);
	}	
}

void basis_dct::plot_basis(ostream & os)
{
	os << "set multiplot layout " << n << "," << n << " rowsfirst scale 1.2;" << endl;

	for(index_t nx = 0; nx < n; nx++)
	for(index_t ny = 0; ny < n; ny++)
	{
		os << "splot v * cos(u), v * sin(u), "; dct(os, nx, ny); os << ";" << endl;
	}
}

void basis_dct::plot_atoms(ostream & os, const a_vec & A)
{
	for(index_t k = 0, nx = 0; nx < n; nx++)
	for(index_t ny = 0; ny < n; ny++, k++)
	{
		os << " + " << A(k) << " * "; dct(os, nx, ny);
	}
}

a_vec basis_dct::dct(const a_vec & x, const a_vec & y, const index_t & nx, const index_t & ny)
{
	return cos( (M_PI * x * nx) / radio ) % cos( (M_PI * y * ny) / radio );
}

void basis_dct::dct(ostream & os, const index_t & nx, const index_t & ny)
{
	os << "cos( (pi * v * cos(u) * " << nx << " ) / " << radio << " ) *";
	os << "cos( (pi * v * sin(u) * " << ny << " ) / " << radio << " )";
}
double basis_dct::get_frequency(size_t idx)
{
	if(idx == 0 ) return INFINITY;
	int tmp = (idx/n > idx%n)? idx/n:idx%n;
	return (2*radio/tmp);
}

} // namespace gproshan::mdict

