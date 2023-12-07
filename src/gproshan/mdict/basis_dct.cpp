#include <gproshan/mdict/basis_dct.h>

#include <cassert>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


basis_dct::basis_dct(const size_t & n, const real_t & r): basis(r, n * n), n_freq(n) {}

void basis_dct::discrete(a_mat & phi, const a_vec & x, const a_vec & y)
{
	assert(phi.n_cols == _dim);

	for(index_t k = 0, nx = 0; nx < n_freq; ++nx)
	for(index_t ny = 0; ny < n_freq; ++ny, ++k)
		phi.col(k) = dct(x, y, nx, ny);
}

void basis_dct::d_discrete(a_mat & phi, const a_vec & x, const a_vec & y, const bool & b)
{
	assert(phi.n_cols == _dim);

	for(index_t k = 0, nx = 0; nx < n_freq; ++nx)
	for(index_t ny = 0; ny < n_freq; ++ny, ++k)
		phi.col(k) = !b ? dct(x, y, nx, ny) : dct(y, x, ny, nx);
}

void basis_dct::plot_basis(std::ostream & os)
{
	os << "set multiplot layout " << n_freq << "," << n_freq << " rowsfirst scale 1.2;" << std::endl;

	for(index_t nx = 0; nx < n_freq; ++nx)
	for(index_t ny = 0; ny < n_freq; ++ny)
	{
		os << "splot v * cos(u), v * sin(u), "; dct(os, nx, ny); os << ";" << std::endl;
	}
}

void basis_dct::plot_atoms(std::ostream & os, const a_vec & A)
{
	for(index_t k = 0, nx = 0; nx < n_freq; ++nx)
	for(index_t ny = 0; ny < n_freq; ++ny, ++k)
	{
		os << " + " << A(k) << " * "; dct(os, nx, ny);
	}
}

a_vec basis_dct::dct(const a_vec & x, const a_vec & y, const index_t nx, const index_t ny)
{
	return cos(M_PI * nx * x / _radio ) % cos(M_PI * ny * y / _radio);
}


a_vec basis_dct::d_dct(const a_vec & x, const a_vec & y, const index_t nx, const index_t ny)
{
	return - (M_PI * nx / _radio) * (sin(M_PI * nx * x / _radio) % cos(M_PI * ny * y / _radio));
}

void basis_dct::dct(std::ostream & os, const index_t nx, const index_t ny)
{
	os << "cos( (pi * v * cos(u) * " << nx << " ) / " << _radio << " ) *";
	os << "cos( (pi * v * sin(u) * " << ny << " ) / " << _radio << " )";
}

real_t basis_dct::freq(const index_t idx)
{
	return !idx ? INFINITY : 2 * _radio / std::max(idx / n_freq, idx % n_freq);
}


} // namespace gproshan::mdict

