#ifndef D_BASIS_COSSINE_H
#define D_BASIS_COSSINE_H

#include "d_basis.h"

class basis_cosine: public basis
{
	private:
		vertex_t r; // rotations

	public:
		basis_cosine(const vertex_t &_radio, const vertex_t & _r);
		void discrete(mat & phi, const mat & xy);
		vec cosine(const vec & x, const vec & y, const vertex_t & c, const vertex_t & alpha);
		void cosine(ostream & os, const vertex_t & c, const vertex_t & alpha);
};

#endif // D_BASIS_COSSINE_H

