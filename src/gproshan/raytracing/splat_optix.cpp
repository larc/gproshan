#include <gproshan/raytracing/splat_optix.h>


//#ifdef GPROSHAN_OPTIX


// geometry processing and shape analysis framework
namespace gproshan::rt {


splat_optix::splat_optix(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats): splat(meshes, model_mats), optix("/src/splat_optix.ptx")
{
	optix_params.traversable = build_as(pointclouds, {mat4::identity()});
	build_sbt();
}


splat_optix::~splat_optix()
{
	cudaFree(d_splats_pcs);
}


} // namespace gproshan

//#endif // GPROSHAN_OPTIX

