#ifdef GPROSHAN_OPTIX

#ifndef RT_OPTIX_H
#define RT_OPTIX_H

#include "che.h"
#include "raytracing.h"

#include <cuda_runtime.h>
#include <optix.h>
#include <optix_stubs.h>

#include <glm/glm.hpp>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


class optix : public raytracing
{
	CUcontext cuda_context;
	CUstream stream;

	OptixDeviceContext optix_context;


	public:
		optix(const std::vector<che *> & meshes);
		~optix();

	private:
		glm::vec4 intersect_li(const glm::vec3 & org, const glm::vec3 & dir, const glm::vec3 & light);
		float intersect_depth(const glm::vec3 & org, const glm::vec3 & dir);
		
		OptixTraversableHandle build_as(const std::vector<che *> & meshes);
		void add_mesh(OptixBuildInput & optix_mesh, uint32_t & optix_trig_flags, const che * mesh);
};


} // namespace gproshan

#endif // RT_OPTIX_H

#endif // GPROSHAN_OPTIX

