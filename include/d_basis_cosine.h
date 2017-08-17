#ifndef D_BASIS_COSSINE_H
#define D_BASIS_COSSINE_H

#include "d_basis.h"

class basis_cosine: public basis
{
	private:
		vertex_t r; // rotations
		vertex_t n; // frequency

	public:
		basis_cosine(const vertex_t &_radio, const size_t & _r, const size_t & _n);
		void discrete(mat & phi, const mat & xy);
	
	private:
		void plot_basis(ostream & os);
		void plot_atoms(ostream & os, const vec & A);
		vec cosine(const vec & x, const vec & y, const vertex_t & c, const vertex_t & alpha);
		void cosine(ostream & os, const vertex_t & c, const vertex_t & alpha);
};

#endif // D_BASIS_COSSINE_H

