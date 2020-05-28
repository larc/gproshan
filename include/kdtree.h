#ifndef KDTREE_H
#define KDTREE_H

#include "vertex.h"


// geometry processing and shape analysis framework
namespace gproshan {


class kdtree
{
	private:
		index_t * nodes = nullptr;

	public:
		kdtree(const vertex * pointcloud, const size_t & n_points);
		virtual ~kdtree();

	private:
		void build(const index_t & n, const vertex * pc, const index_t & i, const index_t & j, const index_t & d);
};


} // namespace gproshan

#endif // KDTREE_H

