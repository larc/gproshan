#ifndef SIMPLIFICATION_H
#define SIMPLIFICATION_H

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


class simplification
{
	private:
		mat4 * Q = nullptr;
		che * mesh = nullptr;
		index_t levels;

	public:
		simplification(che * mesh, const index_t levels_ = 1);
		virtual ~simplification();

		const mat4 & operator [] (const index_t i) const;

	private:
		void execute();
		void compute_quadrics();
		float compute_error(const index_t e);
		void order_edges(index_t * sort_edges, float * error_edges);
		vertex create_vertex(const index_t e);
};


} // namespace gproshan

#endif // SIMPLIFICATION_H

