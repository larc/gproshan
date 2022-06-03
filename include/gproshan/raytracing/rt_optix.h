#ifndef RT_OPTIX_H
#define RT_OPTIX_H

#include <gproshan/mesh/che.h>
#include <gproshan/raytracing/raytracing.h>
#include <gproshan/raytracing/rt_optix_params.h>


#ifdef GPROSHAN_OPTIX

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
	OptixModuleCompileOptions optix_module_compile_opt = {};

	OptixPipeline optix_pipeline;
	OptixPipelineCompileOptions optix_pipeline_compile_opt = {};
	OptixPipelineLinkOptions optix_pipeline_link_opt = {};

	OptixProgramGroup raygen_programs[1];
	OptixProgramGroup miss_programs[2];
	OptixProgramGroup hitgroup_programs[2];

	OptixShaderBindingTable sbt = {};

	launch_params optix_params;
	void * launch_params_buffer = nullptr;

	std::vector<CHE *> dd_mesh;
	std::vector<CHE *> d_mesh;

	void * raygen_records_buffer = nullptr;
	void * miss_records_buffer = nullptr;
	void * hitgroup_records_buffer = nullptr;
	void * as_buffer = nullptr;

	public:
		optix(const std::vector<che *> & meshes, const std::vector<glm::mat4> & model_mats);
		~optix();

		void render(glm::vec4 * img, const render_params & params, const bool & flat);

	private:
		void create_raygen_programs();
		void create_miss_programs();
		void create_hitgroup_programs();
		void create_pipeline();
		void build_sbt();
		OptixTraversableHandle build_as(const std::vector<che *> & meshes, const std::vector<glm::mat4> & model_mats);
		void add_mesh(OptixBuildInput & optix_mesh, CUdeviceptr & d_vertex_ptr, uint32_t & optix_trig_flags, const che * mesh, const glm::mat4 & model_mat);
};


} // namespace gproshan

#endif // GPROSHAN_OPTIX

#endif // RT_OPTIX_H

