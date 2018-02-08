#ifndef D_BASIS_H
#define D_BASIS_H

#include "include.h"

#include <armadillo>
#include <fstream>

using namespace arma;
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
		virtual void discrete(mat & phi, const mat & xy) = 0;
		void plot_basis();
		void plot_atoms(const mat & A);

	private:
		virtual void plot_basis(ostream & os) = 0;
		virtual void plot_atoms(ostream & os, const vec & A) = 0;

	friend class dictionary;
};

} // mdict

#endif // D_BASIS_H

