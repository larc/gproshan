#ifndef DECIMATION_H
#define DECIMATION_H

#include "che.h"

#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

struct corr_t
{
	index_t t;
	vertex alpha;
};

class decimation
{
	private:
		mat * Q;
		che * mesh;
	
	public:
		decimation(che * mesh);
		~decimation();
	
	private:
		void execute();
		void compute_quadrics();
		vertex_t compute_error(const index_t & e);
		void order_edges(index_t * const & sort_edges, vertex_t * const & error_edges);
		vertex create_vertex(const index_t & e);
		corr_t find_corr(const vertex & v, che * mesh, const vector<index_t> & triangles);
};

#endif // DECIMATION_H

