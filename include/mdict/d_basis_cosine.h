#ifndef D_BASIS_COSSINE_H
#define D_BASIS_COSSINE_H

#include "d_basis.h"

// mesh dictionary learning and sparse coding namespace
namespace mdict {

class basis_cosine: public basis
{
	private:
		vertex_t r; // rotations
		vertex_t n; // frequency

	public:
		basis_cosine(const size_t & _r, const size_t & _n, const distance_t & _radio = 0);
		void discrete(a_mat & phi, const a_mat & xy);

	private:
		void plot_basis(ostream & os);
		void plot_atoms(ostream & os, const a_vec & A);
		a_vec cosine(const a_vec & x, const a_vec & y, const vertex_t & c, const vertex_t & alpha);
		void cosine(ostream & os, const vertex_t & c, const vertex_t & alpha);
};

} // mdict

#endif // D_BASIS_COSSINE_H

