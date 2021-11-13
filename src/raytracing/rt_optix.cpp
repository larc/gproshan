#include "raytracing/rt_optix.h"

#ifdef GPROSHAN_OPTIX

#include <cstring>
#include <fstream>
#include <random>

#include <optix_function_table_definition.h>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


void optix_log(unsigned int level, const char * tag, const char * message, void *)
{
	fprintf(stderr, "OptiX [%2u][%12s]: %s\n", level, tag, message);
}

optix::optix(const std::vector<che *> & meshes)
{
	optixInit();

	// create context
	cudaStreamCreate(&stream);

	cuCtxGetCurrent(&cuda_context);

	optixDeviceContextCreate(cuda_context, 0, &optix_context);
	optixDeviceContextSetLogCallback(optix_context, optix_log, nullptr, 4);

	// create module

	OptixModule optix_module;
	OptixModuleCompileOptions optix_module_compile_opt;

	OptixPipeline optix_pipeline;
	OptixPipelineCompileOptions optix_pipeline_compile_opt;
	OptixPipelineLinkOptions optix_pipeline_link_opt;

	optix_module_compile_opt.maxRegisterCount	= 50;
	optix_module_compile_opt.optLevel			= OPTIX_COMPILE_OPTIMIZATION_DEFAULT;
	optix_module_compile_opt.debugLevel			= OPTIX_COMPILE_DEBUG_LEVEL_NONE;

	optix_pipeline_compile_opt							= {};
	optix_pipeline_compile_opt.traversableGraphFlags	= OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_GAS;
	optix_pipeline_compile_opt.usesMotionBlur			= false;
	optix_pipeline_compile_opt.numPayloadValues			= 2;
	optix_pipeline_compile_opt.numAttributeValues		= 2;
	optix_pipeline_compile_opt.exceptionFlags			= OPTIX_EXCEPTION_FLAG_NONE;
	optix_pipeline_compile_opt.pipelineLaunchParamsVariableName = "optix_launch_params";

	optix_pipeline_link_opt.maxTraceDepth		= 2;

	std::ifstream ptx_is(tmp_file_path("rt_optix.ptx"));
	const std::string str_ptx_code = std::string(std::istreambuf_iterator<char>(ptx_is), std::istreambuf_iterator<char>());
	ptx_is.close();

	optixModuleCreateFromPTX(	optix_context,
								&optix_module_compile_opt,
								&optix_pipeline_compile_opt,
								str_ptx_code.c_str(),
								str_ptx_code.size(),
								nullptr, nullptr,	// log message
								&optix_module
								);

	// create programs



	// build as

	build_as(meshes);

	// create pipeline

	// build sbt
}

optix::~optix()
{
}

index_t optix::cast_ray(const glm::vec3 & org, const glm::vec3 & dir)
{
	return NIL;
}

OptixTraversableHandle optix::build_as(const std::vector<che *> & meshes)
{
	OptixTraversableHandle optix_as_handle = {};

	std::vector<OptixBuildInput> optix_meshes(meshes.size());
	std::vector<CUdeviceptr> optix_vertex_ptr(meshes.size());
	std::vector<uint32_t> optix_trig_flags(meshes.size());

	for(index_t i = 0; i < meshes.size(); ++i)
		add_mesh(optix_meshes[i], optix_vertex_ptr[i], optix_trig_flags[i], meshes[i]);

	OptixAccelBuildOptions optix_accel_opt	= {};
	optix_accel_opt.buildFlags 				= OPTIX_BUILD_FLAG_NONE | OPTIX_BUILD_FLAG_ALLOW_COMPACTION;
	optix_accel_opt.motionOptions.numKeys	= 1;
	optix_accel_opt.operation				= OPTIX_BUILD_OPERATION_BUILD;

	OptixAccelBufferSizes optix_gas_buffer_size;
	optixAccelComputeMemoryUsage(	optix_context,
									&optix_accel_opt,
									optix_meshes.data(),
									optix_meshes.size(),
									&optix_gas_buffer_size
									);


	void * d_compacted_size;
	cudaMalloc(&d_compacted_size, sizeof(uint64_t));

	OptixAccelEmitDesc optix_emit_desc;
	optix_emit_desc.type	= OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
	optix_emit_desc.result	= (CUdeviceptr) d_compacted_size;

	void * d_temp_buffer;
	cudaMalloc(&d_temp_buffer, optix_gas_buffer_size.tempSizeInBytes);

	void * d_output_buffer;
	cudaMalloc(&d_output_buffer, optix_gas_buffer_size.outputSizeInBytes);


	optixAccelBuild(	optix_context,
						0,	// stream
						&optix_accel_opt,
						optix_meshes.data(),
						optix_meshes.size(),
						(CUdeviceptr) d_temp_buffer,
						optix_gas_buffer_size.tempSizeInBytes,
						(CUdeviceptr) d_output_buffer,
						optix_gas_buffer_size.outputSizeInBytes,
						&optix_as_handle,
						&optix_emit_desc,
						1
						);

	cudaDeviceSynchronize();

	uint64_t compacted_size;
	cudaMemcpy(&compacted_size, d_compacted_size, sizeof(uint64_t), cudaMemcpyDeviceToHost);

	void * d_as;
	cudaMalloc(&d_as, compacted_size);

	optixAccelCompact(	optix_context,
						0,	// stream
						optix_as_handle,
						(CUdeviceptr) d_as,
						compacted_size,
						&optix_as_handle
						);

	cudaDeviceSynchronize();

	cudaFree(d_output_buffer);
	cudaFree(d_temp_buffer);
	cudaFree(d_compacted_size);

	return optix_as_handle;
}

void optix::add_mesh(OptixBuildInput & optix_mesh, CUdeviceptr & d_vertex_ptr, uint32_t & optix_trig_flags, const che * mesh)
{
	void * d_vertex = nullptr;
	void * d_index = nullptr;

#ifdef GPROSHAN_FLOAT
	cudaMalloc(&d_vertex, mesh->n_vertices * sizeof(vertex));
	cudaMemcpy(d_vertex, &mesh->gt(0), mesh->n_vertices * sizeof(vertex), cudaMemcpyHostToDevice);
#else
	glm::vec3 * vertices = new glm::vec3[mesh->n_vertices];
	cudaMalloc(&d_vertex, mesh->n_vertices * sizeof(float) * 3);

	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices; ++i)
		vertices[i] = glm_vec3(mesh->gt(i));

	cudaMemcpy(d_vertex, vertices, mesh->n_vertices * sizeof(vertex), cudaMemcpyHostToDevice);

	delete [] vertices;
#endif // GPROSHAN_FLOAT

	cudaMalloc(&d_index, mesh->n_half_edges * sizeof(index_t));
	cudaMemcpy(d_index, &mesh->vt(0), mesh->n_half_edges * sizeof(index_t), cudaMemcpyHostToDevice);

	d_vertex_ptr = (CUdeviceptr) d_vertex;

	optix_mesh = {};
	optix_mesh.type = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;

	optix_mesh.triangleArray.vertexFormat			= OPTIX_VERTEX_FORMAT_FLOAT3;
	optix_mesh.triangleArray.vertexStrideInBytes	= 3 * sizeof(float);
	optix_mesh.triangleArray.numVertices			= mesh->n_vertices;
	optix_mesh.triangleArray.vertexBuffers			= &d_vertex_ptr;

	optix_mesh.triangleArray.indexFormat			= OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
	optix_mesh.triangleArray.indexStrideInBytes		= 3 * sizeof(index_t);
	optix_mesh.triangleArray.numIndexTriplets		= mesh->n_faces;
	optix_mesh.triangleArray.indexBuffer			= (CUdeviceptr) d_index;

	optix_trig_flags = 0;

	optix_mesh.triangleArray.flags							= &optix_trig_flags;
	optix_mesh.triangleArray.numSbtRecords					= 1;
	optix_mesh.triangleArray.sbtIndexOffsetBuffer			= 0;
	optix_mesh.triangleArray.sbtIndexOffsetSizeInBytes		= 0;
	optix_mesh.triangleArray.sbtIndexOffsetStrideInBytes	= 0;
}

glm::vec4 optix::intersect_li(const glm::vec3 & org, const glm::vec3 & dir, const glm::vec3 & light, const bool & flat)
{
	return glm::vec4(0.f);
}

float optix::intersect_depth(const glm::vec3 & org, const glm::vec3 & dir)
{
	return 0;
}


} // namespace gproshan

#endif // GPROSHAN_OPTIX

