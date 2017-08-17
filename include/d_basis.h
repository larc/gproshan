#ifndef D_BASIS_H
#define D_BASIS_H

#include "include.h"

#include <armadillo>
#include <fstream>

using namespace arma;
using namespace std;

class basis
{
	protected:
		vertex_t radio;

	public:
		virtual void discrete(mat & phi, const mat & xy) = 0;
		void plot_basis();
		void plot_atoms(const mat & A);
	
	private:
		virtual void plot_basis(ostream & os) = 0;
		virtual void plot_atoms(ostream & os, const vec & A) = 0;
};

#endif // D_BASIS_H

