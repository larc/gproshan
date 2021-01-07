#ifdef GPROSHAN_EMBREE

#ifndef RT_EMBREE_SPLAT_H
#define RT_EMBREE_SPLAT_H

#include "raytracing/rt_embree.h"


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


class embree_splat : public embree
{
	public:
		embree_splat(const std::vector<che *> & meshes, const bool & pointcloud);
	
	private:
		index_t add_pointcloud(const che * mesh);
};


} // namespace gproshan

#endif // RT_EMBREE_SPLAT_H

#endif // GPROSHAN_EMBREE

