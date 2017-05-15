#ifndef DECIMATION_H
#define DECIMATION_H

#include "che.h"

#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

class decimation
{
	private:
		mat * Q;
		che * mesh;
		size_t var_size;
	
	public:
		decimation(che * mesh);
		~decimation();
		void compute_quadrics();
		vertex_t compute_error(const index_t & e);
		void order_edges(index_t * const & sort_edges, vertex_t * const & error_edges);
		vertex create_vertex(const index_t & e);
};


#endif // DECIMATION_H

