#ifndef BASIS_H
#define BASIS_H

#include "include.h"

#include "include_arma.h"
#include <fstream>

using namespace std;


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class dictionary;

class basis
{
	protected:
		distance_t radio;
		size_t dim;

	public:
		virtual ~basis() = default;
		virtual void discrete(a_mat & phi, const a_mat & xy) = 0;
		void plot_basis();
		void plot_atoms(const a_mat & A);
		void plot_patch(const a_mat & A, const a_mat & xyz, index_t i);
		size_t get_dim();
		distance_t get_radio();

	private:
		virtual void plot_basis(std::ostream & os) = 0;
		virtual void plot_atoms(std::ostream & os, const a_vec & A) = 0;

	friend class dictionary;
};


} // namespace gproshan::mdict

#endif // BASIS_H

