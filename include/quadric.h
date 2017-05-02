#ifndef QUADRIC_H
#define QUADRIC_H


#include <string>

#include "include.h"
#include "vertex.h"
#include <armadillo>

using namespace std;
using namespace arma;

class quadric
{
	private:
		mat Q;
		size_t var_size;
	public:
		quadric();
		quadric(vertex normal, vertex p);
		quadric & operator +( quadric & q);
		double eval( vertex v);
	 ~quadric();
};


#endif // CHE_H

