#ifndef D_BASIS_H
#define D_BASIS_H

#include "include.h"

#include "include_arma.h"
#include <fstream>

using namespace std;

// mesh dictionary learning and sparse coding namespace
namespace mdict {

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

	private:
		virtual void plot_basis(std::ostream & os) = 0;
		virtual void plot_atoms(std::ostream & os, const a_vec & A) = 0;

	friend class dictionary;
};

} // mdict

#endif // D_BASIS_H

