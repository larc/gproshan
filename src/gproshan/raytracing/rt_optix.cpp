#include <gproshan/raytracing/rt_optix.h>


#ifdef GPROSHAN_OPTIX


#include <cstring>
#include <fstream>
#include <random>

#include <optix_function_table_definition.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


struct __align__(OPTIX_SBT_RECORD_ALIGNMENT) RaygenRecord
{
	__align__(OPTIX_SBT_RECORD_ALIGNMENT) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
	void * data;
};

struct __align__(OPTIX_SBT_RECORD_ALIGNMENT) MissRecord
{
	__align__(OPTIX_SBT_RECORD_ALIGNMENT) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
	void * data;
};

/*! SBT record for a hitgroup program */
struct __align__(OPTIX_SBT_RECORD_ALIGNMENT) HitgroupRecord
{
	__align__(OPTIX_SBT_RECORD_ALIGNMENT) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
	CHE * data;
};


void optix_log(unsigned int level, const char * tag, const char * message, void *)
{
	fprintf(stderr, "OptiX [%2u][%12s]: %s\n", level, tag, message);
}

optix::optix()
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
	optix_pipeline_compile_opt.pipelineLaunchParamsVariableName = "optix_params";

	optix_pipeline_link_opt.maxTraceDepth = 2;

	std::ifstream ptx_is(std::string(GPROSHAN_DIR) + "/src/rt_optix.ptx");
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

	// create pipeline
	create_pipeline();

	// launch params
	cudaMalloc(&optix_params_buffer, sizeof(launch_params));
}

optix::optix(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats): optix()
{
	// build as
	optix_params.traversable = build_as(meshes, model_mats);

	// build sbt
	build_sbt();
}

optix::~optix()
{
	cudaFree(optix_params_buffer);
	cudaFree(raygen_records_buffer);
	cudaFree(miss_records_buffer);
	cudaFree(hitgroup_records_buffer);
	cudaFree(as_buffer);

	for(index_t i = 0; i < dd_mesh.size(); ++i)
		cuda_free_CHE(dd_mesh[i], d_mesh[i]);
}

void optix::render(vec4 * img, const render_params & params, const bool & flat)
{
	if(params.restart) n_samples = 0;

	optix_params.n_samples = n_samples;
	optix_params.color_buffer = img;

	optix_params.window_width = params.window_width;
	optix_params.window_height = params.window_height;
	if(params.viewport_is_window)
	{
		optix_params.window_width = params.viewport_width;
		optix_params.window_height = params.viewport_height;
	}

	optix_params.viewport_x = params.viewport_x;
	optix_params.viewport_y = params.viewport_y;

	optix_params.flat = flat;
	optix_params.n_lights = params.n_lights;
	memcpy(optix_params.lights, params.lights, sizeof(optix_params.lights));
	memcpy(&optix_params.cam_pos, &params.cam_pos, sizeof(optix_params.cam_pos));
	memcpy(&optix_params.inv_proj_view, &params.inv_proj_view, sizeof(optix_params.inv_proj_view));

	cudaMemcpy(optix_params_buffer, &optix_params, sizeof(launch_params), cudaMemcpyHostToDevice);

	optixLaunch(optix_pipeline,
				stream,
				(CUdeviceptr) optix_params_buffer,
				sizeof(launch_params),
				&sbt,
				params.viewport_width,
				params.viewport_height,
				1
				);

	cudaDeviceSynchronize();

	++n_samples;
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

	if(sizeof_log > 1) gproshan_error_var(log);
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

	if(sizeof_log > 1) gproshan_error_var(log);


	pg_desc.miss.entryFunctionName = "__miss__shadow";

	optixProgramGroupCreate(optix_context,
							&pg_desc,
							1,
							&pg_options,
							log, &sizeof_log,
							&miss_programs[1]
							);

	if(sizeof_log > 1) gproshan_error_var(log);
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

	if(sizeof_log > 1) gproshan_error_var(log);


	pg_desc.hitgroup.entryFunctionNameCH = "__closesthit__shadow";
	pg_desc.hitgroup.entryFunctionNameAH = "__anyhit__shadow";

	optixProgramGroupCreate(optix_context,
							&pg_desc,
							1,
							&pg_options,
							log, &sizeof_log,
							&hitgroup_programs[1]
							);

	if(sizeof_log > 1) gproshan_error_var(log);
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

void optix::build_sbt()
{
	RaygenRecord raygen_records[1];
	for(int i = 0; i < 1; ++i)
	{
		RaygenRecord & rec = raygen_records[i];
		optixSbtRecordPackHeader(raygen_programs[i], &rec);
		rec.data = nullptr;
	}

	cudaMalloc(&raygen_records_buffer, sizeof(RaygenRecord));
	cudaMemcpy(raygen_records_buffer, raygen_records, sizeof(RaygenRecord), cudaMemcpyHostToDevice);
	sbt.raygenRecord = (CUdeviceptr) raygen_records_buffer;


	MissRecord miss_records[2];
	for(int i = 0; i < 2; ++i)
	{
		MissRecord & rec = miss_records[i];
		optixSbtRecordPackHeader(miss_programs[i], &rec);
		rec.data = nullptr;
	}

	cudaMalloc(&miss_records_buffer, 2 * sizeof(MissRecord));
	cudaMemcpy(miss_records_buffer, miss_records, 2 * sizeof(MissRecord), cudaMemcpyHostToDevice);
	sbt.missRecordBase			= (CUdeviceptr) miss_records_buffer;
	sbt.missRecordStrideInBytes	= sizeof(MissRecord);
	sbt.missRecordCount			= 2;


	std::vector<HitgroupRecord> hitgroup_records;
	for(index_t i = 0; i < d_mesh.size(); ++i)
	for(index_t r = 0; r < 2; ++r)
	{
		HitgroupRecord rec;
		optixSbtRecordPackHeader(hitgroup_programs[r], &rec);
		rec.data = d_mesh[i];
		hitgroup_records.push_back(rec);
	}

	cudaMalloc(&hitgroup_records_buffer, hitgroup_records.size() * sizeof(HitgroupRecord));
	cudaMemcpy(hitgroup_records_buffer, hitgroup_records.data(), hitgroup_records.size() * sizeof(HitgroupRecord), cudaMemcpyHostToDevice);
	sbt.hitgroupRecordBase			= (CUdeviceptr) hitgroup_records_buffer;
	sbt.hitgroupRecordStrideInBytes	= sizeof(HitgroupRecord);
	sbt.hitgroupRecordCount			= hitgroup_records.size();
}

OptixTraversableHandle optix::build_as(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats)
{
	OptixTraversableHandle optix_as_handle = {};

	std::vector<OptixBuildInput> optix_meshes(meshes.size());
	std::vector<CUdeviceptr> optix_vertex_ptr(meshes.size());
	std::vector<uint32_t> optix_trig_flags(meshes.size());

	for(index_t i = 0; i < meshes.size(); ++i)
		add_mesh(optix_meshes[i], optix_vertex_ptr[i], optix_trig_flags[i], meshes[i], model_mats[i]);

	OptixAccelBuildOptions optix_accel_opt	=	{};
	optix_accel_opt.buildFlags 				=	OPTIX_BUILD_FLAG_ALLOW_RANDOM_VERTEX_ACCESS |
												OPTIX_BUILD_FLAG_ALLOW_COMPACTION;
	optix_accel_opt.operation				=	OPTIX_BUILD_OPERATION_BUILD;

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

	cudaMalloc(&as_buffer, compacted_size);

	optixAccelCompact(	optix_context,
						0,	// stream
						optix_as_handle,
						(CUdeviceptr) as_buffer,
						compacted_size,
						&optix_as_handle
						);

	cudaDeviceSynchronize();

	cudaFree(d_output_buffer);
	cudaFree(d_temp_buffer);
	cudaFree(d_compacted_size);

	return optix_as_handle;
}

void optix::add_mesh(OptixBuildInput & optix_mesh, CUdeviceptr & d_vertex_ptr, uint32_t & optix_trig_flags, const che * mesh, const mat4 & model_mat)
{
	CHE * dd_m, * d_m;
	CHE h_m(mesh);

	if(mesh->is_scene())
	{
		h_m.n_half_edges = mesh->n_vertices;
		h_m.n_trigs = mesh->n_vertices / 3;
		h_m.VT = new index_t[mesh->n_vertices];

		#pragma omp parallel for
		for(index_t i = 0; i < mesh->n_vertices; ++i)
			h_m.VT[i] = i;
	}

	cuda_create_CHE(&h_m, dd_m, d_m, true, true);
	dd_mesh.push_back(dd_m);
	d_mesh.push_back(d_m);

	float * d_model_mat = nullptr;
	cudaMalloc(&d_model_mat, sizeof(model_mat));
	cudaMemcpy(d_model_mat, &model_mat, sizeof(model_mat), cudaMemcpyHostToDevice);

	d_vertex_ptr = (CUdeviceptr) dd_m->GT;

	optix_mesh = {};
	optix_mesh.type = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;

	optix_mesh.triangleArray.vertexFormat			= OPTIX_VERTEX_FORMAT_FLOAT3;
	optix_mesh.triangleArray.vertexStrideInBytes	= 3 * sizeof(float);
	optix_mesh.triangleArray.numVertices			= h_m.n_vertices;
	optix_mesh.triangleArray.vertexBuffers			= &d_vertex_ptr;

	optix_mesh.triangleArray.indexFormat			= OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
	optix_mesh.triangleArray.indexStrideInBytes		= 3 * sizeof(index_t);
	optix_mesh.triangleArray.numIndexTriplets		= h_m.n_trigs;
	optix_mesh.triangleArray.indexBuffer			= (CUdeviceptr) dd_m->VT;

	optix_mesh.triangleArray.transformFormat		= OPTIX_TRANSFORM_FORMAT_MATRIX_FLOAT12;
	optix_mesh.triangleArray.preTransform			= (CUdeviceptr) d_model_mat;

	optix_trig_flags = 0;

	optix_mesh.triangleArray.flags							= &optix_trig_flags;
	optix_mesh.triangleArray.numSbtRecords					= 1;
	optix_mesh.triangleArray.sbtIndexOffsetBuffer			= 0;
	optix_mesh.triangleArray.sbtIndexOffsetSizeInBytes		= 0;
	optix_mesh.triangleArray.sbtIndexOffsetStrideInBytes	= 0;

	if(mesh->is_scene())
		delete [] h_m.VT;
}


} // namespace gproshan

#endif // GPROSHAN_OPTIX

