#ifdef GPROSHAN_OPTIX

#ifndef RT_OPTIX_H
#define RT_OPTIX_H

#include "mesh/che.h"
#include "raytracing/raytracing.h"

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
	OptixModule optix_module;
	OptixProgramGroup raygen_programs[1];
	OptixProgramGroup miss_programs[2];
	OptixProgramGroup hitgroup_programs[2];

	public:
		optix(const std::vector<che *> & meshes);
		~optix();

		virtual index_t cast_ray(const glm::vec3 & org, const glm::vec3 & dir);

	private:
		glm::vec4 intersect_li(const glm::vec3 & org, const glm::vec3 & dir, const glm::vec3 & light, const bool & flat);
		float intersect_depth(const glm::vec3 & org, const glm::vec3 & dir);

		void create_raygen_programs();
		void create_miss_programs();
		void create_hitgroup_programs();
		OptixTraversableHandle build_as(const std::vector<che *> & meshes);
		void add_mesh(OptixBuildInput & optix_mesh, CUdeviceptr & d_vertex_ptr, uint32_t & optix_trig_flags, const che * mesh);
};


} // namespace gproshan

#endif // RT_OPTIX_H

#endif // GPROSHAN_OPTIX

