#ifndef DECIMATION_H
#define DECIMATION_H

#include "che.h"

#include <string>
#include "include_arma.h"

using namespace std;

class decimation
{
	private:
		a_mat * Q;
		che * mesh;
		corr_t * corr;
		index_t levels;

	public:
		decimation(che * mesh, const vertex *const & normals, const index_t & levels_ = 1);
		~decimation();
		operator const corr_t * ();

	private:
		void execute(const vertex *const & normals);
		void compute_quadrics();
		vertex_t compute_error(const index_t & e);
		void order_edges(index_t * const & sort_edges, vertex_t * const & error_edges);
		vertex create_vertex(const index_t & e);
};

#endif // DECIMATION_H

