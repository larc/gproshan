#ifndef D_BASIS_DCT_H
#define D_BASIS_DCT_H

#include "d_basis.h"

// mesh dictionary learning and sparse coding namespace
namespace mdict {

class basis_dct: public basis
{
	private:
		vertex_t n; // frequency

	public:
		basis_dct(const size_t & _n, const distance_t & _radio = 0);
		void discrete(a_mat & phi, const a_mat & xy);

	private:
		void plot_basis(ostream & os);
		void plot_atoms(ostream & os, const a_vec & A);
		a_vec dct(const a_vec & x, const a_vec & y, const index_t & nx, const index_t & ny);
		void dct(ostream & os, const index_t & nx, const index_t & ny);
};

} // mdict

#endif // D_BASIS_DCT_H

