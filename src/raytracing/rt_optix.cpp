#ifdef GPROSHAN_OPTIX

#include "raytracing/rt_optix.h"

#include "mesh/che.cuh"

#include <cstring>
#include <fstream>
#include <random>

#include <optix_function_table_definition.h>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


/*! SBT record for a raygen program */
struct __align__(OPTIX_SBT_RECORD_ALIGNMENT) RaygenRecord
{
	__align__(OPTIX_SBT_RECORD_ALIGNMENT) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
	// just a dummy value - later examples will use more interesting
	// data here
	void * data;
};

/*! SBT record for a miss program */
struct __align__(OPTIX_SBT_RECORD_ALIGNMENT) MissRecord
{
	__align__(OPTIX_SBT_RECORD_ALIGNMENT) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
	// just a dummy value - later examples will use more interesting
	// data here
	void * data;
};

/*! SBT record for a hitgroup program */
struct __align__(OPTIX_SBT_RECORD_ALIGNMENT) HitgroupRecord
{
	__align__(OPTIX_SBT_RECORD_ALIGNMENT) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
	mesh_sbt_data mesh;
};


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

	optix_module_compile_opt.maxRegisterCount	= 50;
	optix_module_compile_opt.optLevel			= OPTIX_COMPILE_OPTIMIZATION_DEFAULT;
	optix_module_compile_opt.debugLevel			= OPTIX_COMPILE_DEBUG_LEVEL_NONE;

	optix_pipeline_compile_opt							= {};
	optix_pipeline_compile_opt.traversableGraphFlags	= OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_GAS;
	optix_pipeline_compile_opt.usesMotionBlur			= false;
	optix_pipeline_compile_opt.numPayloadValues			= 2;
	optix_pipeline_compile_opt.numAttributeValues		= 2;
	optix_pipeline_compile_opt.exceptionFlags			= OPTIX_EXCEPTION_FLAG_NONE;
	optix_pipeline_compile_opt.pipelineLaunchParamsVariableName = "params";

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
	create_raygen_programs();
	create_miss_programs();
	create_hitgroup_programs();

	// build as
	build_as(meshes);

	// create pipeline
	create_pipeline();
	
	// build sbt


}

optix::~optix()
{
	for(index_t i = 0; i < dd_mesh.size(); ++i)
		cuda_free_CHE(dd_mesh[i], d_mesh[i]);
}

index_t optix::cast_ray(const glm::vec3 & org, const glm::vec3 & dir)
{
	return NIL;
}

void optix::create_raygen_programs()
{
	char log[2048];
	size_t sizeof_log = sizeof(log);

	OptixProgramGroupOptions pg_options	= {};
	OptixProgramGroupDesc pg_desc		= {};
	pg_desc.kind						= OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
	pg_desc.raygen.module				= optix_module;
	pg_desc.raygen.entryFunctionName	= "__raygen__render_frame";

	optixProgramGroupCreate(optix_context,
							&pg_desc,
							1,
							&pg_options,
							log, &sizeof_log,
							&raygen_programs[0]
							);

	if(sizeof_log > 1) gproshan_log_var(log);
}

void optix::create_miss_programs()
{
	char log[2048];
	size_t sizeof_log = sizeof(log);

	OptixProgramGroupOptions pg_options	= {};
	OptixProgramGroupDesc pg_desc		= {};
	pg_desc.kind						= OPTIX_PROGRAM_GROUP_KIND_MISS;
	pg_desc.miss.module					= optix_module;


	pg_desc.miss.entryFunctionName = "__miss__radiance";

	optixProgramGroupCreate(optix_context,
							&pg_desc,
							1,
							&pg_options,
							log, &sizeof_log,
							&miss_programs[0]
							);

	if(sizeof_log > 1) gproshan_log_var(log);


	pg_desc.miss.entryFunctionName = "__miss__shadow";

	optixProgramGroupCreate(optix_context,
							&pg_desc,
							1,
							&pg_options,
							log, &sizeof_log,
							&miss_programs[1]
							);

	if(sizeof_log > 1) gproshan_log_var(log);
}

void optix::create_hitgroup_programs()
{
	char log[2048];
	size_t sizeof_log = sizeof(log);

	OptixProgramGroupOptions pg_options	= {};
	OptixProgramGroupDesc pg_desc		= {};
	pg_desc.kind						= OPTIX_PROGRAM_GROUP_KIND_HITGROUP;
	pg_desc.hitgroup.moduleCH			= optix_module;
	pg_desc.hitgroup.moduleAH			= optix_module;


	pg_desc.hitgroup.entryFunctionNameCH = "__closesthit__radiance";
	pg_desc.hitgroup.entryFunctionNameAH = "__anyhit__radiance";

	optixProgramGroupCreate(optix_context,
							&pg_desc,
							1,
							&pg_options,
							log, &sizeof_log,
							&hitgroup_programs[0]
							);

	if(sizeof_log > 1) gproshan_log_var(log);


	pg_desc.hitgroup.entryFunctionNameCH = "__closesthit__shadow";
	pg_desc.hitgroup.entryFunctionNameAH = "__anyhit__shadow";

	optixProgramGroupCreate(optix_context,
							&pg_desc,
							1,
							&pg_options,
							log, &sizeof_log,
							&hitgroup_programs[1]
							);

	if(sizeof_log > 1) gproshan_log_var(log);
}

void optix::create_pipeline()
{
	std::vector<OptixProgramGroup> program_groups;
	program_groups.push_back(raygen_programs[0]);
	program_groups.push_back(hitgroup_programs[0]);
	program_groups.push_back(hitgroup_programs[1]);
	program_groups.push_back(miss_programs[0]);
	program_groups.push_back(miss_programs[1]);

	char log[2048];
	size_t sizeof_log = sizeof(log);
	
	optixPipelineCreate(optix_context,
						&optix_pipeline_compile_opt,
						&optix_pipeline_link_opt,
						program_groups.data(),
						program_groups.size(),
						log, &sizeof_log,
						&optix_pipeline
						);
	
	if(sizeof_log > 1) gproshan_log_var(log);

	optixPipelineSetStackSize(optix_pipeline, 2 * 1024, 2 * 1024, 2 * 1024, 1);
	
	if(sizeof_log > 1) gproshan_log_var(log);
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
	CHE * dd_m, * d_m;
	CHE h_m(mesh);

	cuda_create_CHE(&h_m, dd_m, d_m);
	dd_mesh.push_back(dd_m);
	d_mesh.push_back(d_m);

	
	d_vertex_ptr = (CUdeviceptr) dd_m->GT;

	optix_mesh = {};
	optix_mesh.type = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;

	optix_mesh.triangleArray.vertexFormat			= OPTIX_VERTEX_FORMAT_FLOAT3;
	optix_mesh.triangleArray.vertexStrideInBytes	= sizeof(vertex);
	optix_mesh.triangleArray.numVertices			= mesh->n_vertices;
	optix_mesh.triangleArray.vertexBuffers			= &d_vertex_ptr;

	optix_mesh.triangleArray.indexFormat			= OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
	optix_mesh.triangleArray.indexStrideInBytes		= 3 * sizeof(index_t);
	optix_mesh.triangleArray.numIndexTriplets		= mesh->n_faces;
	optix_mesh.triangleArray.indexBuffer			= (CUdeviceptr) dd_m->VT;

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

