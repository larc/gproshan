#ifndef DECIMATION_H
#define DECIMATION_H

#include "che.h"
#include "include_arma.h"

#include <string>


// geometry processing and shape analysis framework
namespace gproshan {


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
		real_t compute_error(const index_t & e);
		void order_edges(index_t * const & sort_edges, real_t * const & error_edges);
		vertex create_vertex(const index_t & e);
};


} // namespace gproshan

#endif // DECIMATION_H

