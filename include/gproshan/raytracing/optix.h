#ifndef RT_OPTIX_H
#define RT_OPTIX_H

#include <gproshan/mesh/che_cuda.h>
#include <gproshan/scenes/scene.h>
#include <gproshan/raytracing/raytracing.h>
#include <gproshan/raytracing/optix_params.h>


#ifdef GPROSHAN_OPTIX

#include <cuda_runtime.h>
#include <optix.h>
#include <optix_stubs.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


class optix : public raytracing
{
	protected:
		CUcontext cuda_context;
		CUstream stream;

		OptixDeviceContext _context;

		OptixModule _module;
		OptixModuleCompileOptions _module_compile_opt = {};

		OptixPipeline _pipeline;
		OptixPipelineCompileOptions _pipeline_compile_opt = {};
		OptixPipelineLinkOptions _pipeline_link_opt = {};

		OptixProgramGroup raygen_programs[1];
		OptixProgramGroup miss_programs[2];
		OptixProgramGroup hitgroup_programs[2];

		OptixShaderBindingTable sbt = {};

		optix_params params;
		optix_params * params_buffer = nullptr;

		std::vector<che *> d_mesh;

		void * raygen_records_buffer = nullptr;
		void * miss_records_buffer = nullptr;
		void * hitgroup_records_buffer = nullptr;
		void * as_buffer = nullptr;

		std::vector<unsigned char *> tex_data;

	public:
		optix(const std::string & program = "optix.optixir");
		optix(const std::vector<const che *> & meshes, const std::vector<mat4> & model_mats);
		virtual ~optix();

		virtual void render(vec4 * img, const render_params & params, const bool flat);

	protected:
		void create_raygen_programs();
		void create_miss_programs();
		void create_hitgroup_programs();
		void create_pipeline();
		void build_sbt();
		OptixTraversableHandle build_as(const std::vector<const che *> & meshes, const std::vector<mat4> & model_mats);
		void add_mesh(OptixBuildInput & _mesh, CUdeviceptr & d_vertex_ptr, uint32_t & _trig_flags, const che * mesh, const mat4 & model_mat);
};


} // namespace gproshan

#endif // GPROSHAN_OPTIX

#endif // RT_OPTIX_H

