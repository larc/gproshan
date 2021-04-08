#ifndef SIMPLIFICATION_H
#define SIMPLIFICATION_H

#include "mesh/che.h"
#include "include_arma.h"

#include <string>


// geometry processing and shape analysis framework
namespace gproshan {


class simplification
{
	private:
		a_mat * Q;
		che * mesh;
		index_t levels;

	public:
		simplification(che * mesh, const index_t & levels_ = 1);
		~simplification();

	private:
		void execute();
		void compute_quadrics();
		real_t compute_error(const index_t & e);
		void order_edges(index_t * const & sort_edges, real_t * const & error_edges);
		vertex create_vertex(const index_t & e);
};


} // namespace gproshan

#endif // SIMPLIFICATION_H

