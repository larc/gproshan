#ifndef BASIS_H
#define BASIS_H

#include <gproshan/include.h>
#include <gproshan/include_arma.h>

#include <fstream>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class basis
{
	protected:
		real_t _radio;
		size_t _dim;

	public:
		basis(const real_t r, const size_t d);
		virtual ~basis() = default;
		virtual void discrete(a_mat & phi, const a_vec & x, const a_vec & y) = 0;
		virtual void d_discrete(a_mat & phi, const a_vec & x, const a_vec & y, const bool & b) = 0;
		virtual real_t freq(const index_t idx) = 0;
		real_t & radio();
		size_t dim() const;
		void plot_basis();
		void plot_atoms(const a_mat & A);
		void plot_patch(const a_mat & A, const a_mat & xyz, const index_t p);

	private:
		virtual void plot_basis(std::ostream & os) = 0;
		virtual void plot_atoms(std::ostream & os, const a_vec & A) = 0;
};


} // namespace gproshan::mdict

#endif // BASIS_H

