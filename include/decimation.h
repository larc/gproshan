#ifndef DECIMATION_H
#define DECIMATION_H

#include "che.h"

#include <string>
#include <armadillo>
#include <queue>

using namespace std;
using namespace arma;

class err_edge
{
	public:
	vertex_t error;
	index_t edge;
};

bool operator > ( err_edge & a, err_edge & b )
{
		return a.error < b.error;
}
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
		priority_queue< err_edge > order_edges();
		vertex  create_vertex(const index_t & e);
};


#endif // DECIMATION_H

