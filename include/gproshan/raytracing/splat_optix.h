#ifndef RT_SPLAT_OPTIX_H
#define RT_SPLAT_OPTIX_H

#ifdef GPROSHAN_OPTIX

#include <gproshan/raytracing/splat.h>
#include <gproshan/raytracing/optix.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


class splat_optix: public splat, public optix
{
	public:
		splat_optix(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats): splat(meshes, model_mats), optix("/src/splat_optix.ptx")
		{
			optix_params.traversable = build_as(pointclouds, {mat4::identity()});
			build_sbt();
		}
};


} // namespace gproshan

#endif // GPROSHAN_OPTIX

#endif // RT_SPLAT_EMBREE_H

