#ifndef D_BASIS_COSSINE_H
#define D_BASIS_COSSINE_H

#include "d_basis.h"

// mesh dictionary learning and sparse coding namespace
namespace mdict {

class basis_cosine: public basis
{
	private:
		real_t r; // rotations
		real_t n; // frequency

	public:
		basis_cosine(const size_t & _r, const size_t & _n, const distance_t & _radio = 0);
		void discrete(a_mat & phi, const a_mat & xy);

	private:
		void plot_basis(std::ostream & os);
		void plot_atoms(std::ostream & os, const a_vec & A);
		a_vec cosine(const a_vec & x, const a_vec & y, const real_t & c, const real_t & alpha);
		void cosine(std::ostream & os, const real_t & c, const real_t & alpha);
};

} // mdict

#endif // D_BASIS_COSSINE_H

