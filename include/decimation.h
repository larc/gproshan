#ifndef DECIMATION_H
#define DECIMATION_H


#include <string>

#include "include.h"
#include "vertex.h"
#include <armadillo>
#include "che.h"


using namespace std;
using namespace arma;

class decimation
{
	private:
		mat *Q;
		size_t var_size;
	public:
		decimation(che * mesh, size_t n_choosed_edges);
		double eval( vertex_t v);
		double compute_error(const vertex_t & a, const vertex_t & b);
	 ~decimation();
};


#endif // CHE_H

