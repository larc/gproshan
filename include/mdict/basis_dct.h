#ifndef BASIS_DCT_H
#define BASIS_DCT_H

#include "basis.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class basis_dct: public basis
{
	private:
		real_t n; // frequency

	public:
		basis_dct(const size_t & _n, const distance_t & _radio = 0);
		void discrete(a_mat & phi, const a_mat & xy);

	private:
		void plot_basis(std::ostream & os);
		void plot_atoms(std::ostream & os, const a_vec & A);
		a_vec dct(const a_vec & x, const a_vec & y, const index_t & nx, const index_t & ny);
		void dct(std::ostream & os, const index_t & nx, const index_t & ny);
};


} // namespace gproshan::mdict

#endif // BASIS_DCT_H

