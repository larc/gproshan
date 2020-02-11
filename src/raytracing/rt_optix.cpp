#include "rt_optix.h"

#ifdef GPROSHAN_OPTIX

#include <iostream>
#include <random>
#include <cstring>

#include <optix_function_table_definition.h>

// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


optix::optix()
{
	optixInit();
}

optix::~optix()
{
}

void optix::build_bvh()
{
}

index_t optix::add_sphere(const glm::vec4 & xyzr)
{
	return 0;
}

index_t optix::add_mesh(const che * mesh, const glm::mat4 & model_matrix)
{
	return 0;
}

index_t optix::add_point_cloud(const che * mesh, const glm::mat4 & model_matrix)
{
	return 0;
}

const glm::vec4 optix::intersect_li(const glm::vec3 & org, const glm::vec3 & dir, const glm::vec3 & light)
{
	return glm::vec4(0.f);
}

const float optix::intersect_depth(const glm::vec3 & org, const glm::vec3 & dir)
{
	return 0;
}


} // namespace gproshan

#endif // GPROSHAN_OPTIX

