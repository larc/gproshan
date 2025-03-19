#ifndef RT_SPLAT_EMBREE_H
#define RT_SPLAT_EMBREE_H

#include <gproshan/raytracing/splat.h>
#include <gproshan/raytracing/embree.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


class splat_embree: public splat, public embree
{
	public:
		splat_embree(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats);

		bool closesthit_radiance(	vertex & color,
									vertex & attenuation,
									vertex & position,
									vertex & ray_dir,
									float & dist,
									random<float> & rnd,
									const render_params & params,
									const bool
									) const;
};


} // namespace gproshan

#endif // RT_SPLAT_EMBREE_H

