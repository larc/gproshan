#ifndef DECIMATION_H
#define DECIMATION_H


#include <string>

#include "include.h"
#include "vertex.h"
#include <armadillo>

using namespace std;
using namespace arma;

class decimation
{
	private:
		mat Q;
		size_t var_size;
	public:
		decimation();
		decimation(vertex normal, vertex p);
		double eval( vertex v);
		double compute_error(const vertex_t & a, const vertex_t & b);
	 ~decimation();
};


#endif // CHE_H

