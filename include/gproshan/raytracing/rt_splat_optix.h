#ifndef RT_SPLAT_OPTIX_H
#define RT_SPLAT_OPTIX_H

#include <gproshan/raytracing/splat.h>
#include <gproshan/raytracing/rt_optix.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


class splat_optix: public splat, public optix
{
	public:
		splat_optix(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats): splat(meshes, model_mats)
		{
			// build as
			//optix_params.traversable = build_as(meshes, {mat4::identity()});

			// build sbt
			//build_sbt();
		}
};


} // namespace gproshan

#endif // RT_SPLAT_EMBREE_H

