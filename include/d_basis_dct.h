#ifndef D_BASIS_DCT_H
#define D_BASIS_DCT_H

#include "d_basis.h"

class basis_dct: public basis
{
	private:
		vertex_t n; // frequency

	public:
		basis_dct(const vertex_t &_radio, const size_t & _n);
		void discrete(mat & phi, const mat & xy);
	
	private:
		void plot_basis(ostream & os);
		void plot_atoms(ostream & os, const vec & A);
		vec dct(const vec & x, const vec & y, const index_t & nx, const index_t & ny);
		void dct(ostream & os, const index_t & nx, const index_t & ny);
};

#endif // D_BASIS_DCT_H

