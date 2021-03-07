#include "mesh/kdtree.h"

#include <algorithm>


// geometry processing and shape analysis framework
namespace gproshan {


kdtree::kdtree(const vertex * pointcloud, const size_t & n_points)
{
	nodes = new index_t[n_points >> 1];

	build(0, pointcloud, 0, n_points, 0);
}

kdtree::~kdtree()
{
	delete [] nodes;
}

void kdtree::build(const index_t & n, const vertex * pc, const index_t & i, const index_t & j, const index_t & d)
{

	if(i == j)
	{
		nodes[i] = i;
		return;
	}
}


} // namespace gproshan

