#ifndef SIMPLIFICATION_H
#define SIMPLIFICATION_H

#include <gproshan/mesh/che.h>

#include <string>
#include <armadillo>


// geometry processing and shape analysis framework
namespace gproshan {


class simplification
{
	private:
		arma::fmat * Q;
		che * mesh;
		index_t levels;

	public:
		simplification(che * mesh, const index_t levels_ = 1);
		~simplification();

	private:
		void execute();
		void compute_quadrics();
		real_t compute_error(const index_t e);
		void order_edges(index_t * sort_edges, real_t * error_edges);
		vertex create_vertex(const index_t e);
};


} // namespace gproshan

#endif // SIMPLIFICATION_H

