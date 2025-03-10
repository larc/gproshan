#include <gproshan/raytracing/optix.h>


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
	che * data;
};


void optix_log(index_t level, const char * tag, const char * message, void *)
{
	fprintf(stderr, "OptiX [%2u][%12s]: %s\n", level, tag, message);
}

optix::optix(const std::string & program, const unsigned int nthreads)
{
	optixInit();

	cudaStreamCreate(&stream);

	cuCtxGetCurrent(&cuda_context);

	optixDeviceContextCreate(cuda_context, 0, &_context);
	optixDeviceContextSetLogCallback(_context, optix_log, nullptr, 4);

	_module_compile_opt.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_NONE;

	_pipeline_compile_opt.traversableGraphFlags		= OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_GAS;
	_pipeline_compile_opt.usesMotionBlur			= false;
	_pipeline_compile_opt.numPayloadValues			= 4;
	_pipeline_compile_opt.numAttributeValues		= 4;
	_pipeline_compile_opt.exceptionFlags			= OPTIX_EXCEPTION_FLAG_NONE;
	_pipeline_compile_opt.pipelineLaunchParamsVariableName = "params";

	_pipeline_link_opt.maxTraceDepth = 2;

	std::ifstream is(tmp_file_path(program));
	const std::string program_src = std::string(std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>());
	is.close();

	optixModuleCreate(	_context,
						&_module_compile_opt,
						&_pipeline_compile_opt,
						program_src.c_str(),
						size(program_src),
						nullptr, nullptr,	// log message
						&_module
						);

	create_raygen_programs();
	create_miss_programs();
	create_hitgroup_programs();

	create_pipeline();

	params_buffer.assign(nthreads, nullptr);
	gproshan_error_var(size(params_buffer));
	for(auto & p: params_buffer)
		cudaMalloc(&p, sizeof(optix_params));
}

optix::optix(const std::vector<const che *> & meshes, const std::vector<mat4> & model_mats): optix()
{
	params.traversable = build_as(meshes, model_mats);
	build_sbt();
}

optix::~optix()
{
	for(auto & p: params_buffer)
		cudaFree(p);

	cudaFree(raygen_records_buffer);
	cudaFree(miss_records_buffer);
	cudaFree(hitgroup_records_buffer);
	cudaFree(as_buffer);

	for(index_t i = 0; i < size(d_mesh); ++i)
		delete d_mesh[i];

	cudaFree(params.sc.materials);
	cudaFree(params.sc.textures);
	cudaFree(params.sc.trig_mat);
	cudaFree(params.sc.texcoords);

	for(unsigned char * data: tex_data)
		cudaFree(data);
}

void optix::render(vec4 * img, const render_params & rp, const bool flat)
{
	update_params(img, rp, flat);
	render(rp.thread, rp.viewport_size.x() * rp.viewport_size.y());
}

void optix::update_params(vec4 * img, const render_params & rp, const bool flat)
{
	optix_params tmp_params = params;

	tmp_params.depth = rp.depth;
	tmp_params.thread = rp.thread;
	tmp_params.n_frames = rp.n_frames;
	tmp_params.n_samples = rp.n_samples;
	tmp_params.color_buffer = img;

	tmp_params.viewport_size = rp.viewport_size;
	tmp_params.window_size = rp.window_size;
	if(rp.viewport_is_window)
		tmp_params.window_size = rp.viewport_size;

	tmp_params.viewport_pos = rp.viewport_pos;

	tmp_params.flat = flat;
	tmp_params.cam_pos = rp.cam_pos;
	tmp_params.inv_proj_view = rp.inv_proj_view;
	tmp_params.ambient = rp.ambient;
	tmp_params.n_lights = rp.n_lights;
	memcpy(tmp_params.lights, rp.lights, sizeof(params.lights));

	cudaMemcpy(params_buffer[rp.thread], &tmp_params, sizeof(optix_params), cudaMemcpyHostToDevice);
}

void optix::render(const unsigned int thread, const unsigned int nrays)
{
	optixLaunch(_pipeline
				, stream
				, (CUdeviceptr) params_buffer[thread]
				, sizeof(optix_params)
				, &sbt
				, nrays
				, 1
				, 1
				);

	cudaDeviceSynchronize();
}


void optix::create_raygen_programs()
{
	char log[2048];
	size_t sizeof_log = sizeof(log);

	OptixProgramGroupOptions pg_options	= {};
	OptixProgramGroupDesc pg_desc		= {};
	pg_desc.kind						= OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
	pg_desc.raygen.module				= _module;
	pg_desc.raygen.entryFunctionName	= "__raygen__render_frame";

	optixProgramGroupCreate(_context,
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
	pg_desc.miss.module					= _module;


	pg_desc.miss.entryFunctionName = "__miss__radiance";

	optixProgramGroupCreate(_context,
							&pg_desc,
							1,
							&pg_options,
							log, &sizeof_log,
							&miss_programs[0]
							);

	if(sizeof_log > 1) gproshan_error_var(log);


	pg_desc.miss.entryFunctionName = "__miss__shadow";

	optixProgramGroupCreate(_context,
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
	pg_desc.hitgroup.moduleCH			= _module;
	pg_desc.hitgroup.moduleAH			= _module;


	pg_desc.hitgroup.entryFunctionNameCH = "__closesthit__radiance";
	pg_desc.hitgroup.entryFunctionNameAH = "__anyhit__radiance";

	optixProgramGroupCreate(_context,
							&pg_desc,
							1,
							&pg_options,
							log, &sizeof_log,
							&hitgroup_programs[0]
							);

	if(sizeof_log > 1) gproshan_error_var(log);


	pg_desc.hitgroup.entryFunctionNameCH = "__closesthit__shadow";
	pg_desc.hitgroup.entryFunctionNameAH = "__anyhit__shadow";

	optixProgramGroupCreate(_context,
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

	optixPipelineCreate(_context,
						&_pipeline_compile_opt,
						&_pipeline_link_opt,
						program_groups.data(),
						size(program_groups),
						log, &sizeof_log,
						&_pipeline
						);

	optixPipelineSetStackSize(_pipeline, 2 * 1024, 2 * 1024, 2 * 1024, 1);
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
	for(index_t i = 0; i < size(d_mesh); ++i)
	for(index_t r = 0; r < 2; ++r)
	{
		HitgroupRecord rec;
		optixSbtRecordPackHeader(hitgroup_programs[r], &rec);
		che_cuda & m = *(che_cuda *) d_mesh[i];
		rec.data = m;
		hitgroup_records.push_back(rec);
	}

	cudaMalloc(&hitgroup_records_buffer, size(hitgroup_records) * sizeof(HitgroupRecord));
	cudaMemcpy(hitgroup_records_buffer, hitgroup_records.data(), size(hitgroup_records) * sizeof(HitgroupRecord), cudaMemcpyHostToDevice);
	sbt.hitgroupRecordBase			= (CUdeviceptr) hitgroup_records_buffer;
	sbt.hitgroupRecordStrideInBytes	= sizeof(HitgroupRecord);
	sbt.hitgroupRecordCount			= size(hitgroup_records);
}

OptixTraversableHandle optix::build_as(const std::vector<const che *> & meshes, const std::vector<mat4> & model_mats)
{
	OptixTraversableHandle _as_handle = {};

	std::vector<OptixBuildInput> _meshes(size(meshes));
	std::vector<CUdeviceptr> _vertex_ptr(size(meshes));
	std::vector<uint32_t> _trig_flags(size(meshes));

	for(index_t i = 0; i < size(meshes); ++i)
		add_mesh(_meshes[i], _vertex_ptr[i], _trig_flags[i], meshes[i], model_mats[i]);

	OptixAccelBuildOptions _accel_opt	= {};
	_accel_opt.buildFlags 				= OPTIX_BUILD_FLAG_ALLOW_RANDOM_VERTEX_ACCESS |
											OPTIX_BUILD_FLAG_ALLOW_COMPACTION;
	_accel_opt.operation				= OPTIX_BUILD_OPERATION_BUILD;

	OptixAccelBufferSizes _gas_buffer_size;
	optixAccelComputeMemoryUsage(	_context,
									&_accel_opt,
									_meshes.data(),
									size(_meshes),
									&_gas_buffer_size
									);


	void * d_compacted_size;
	cudaMalloc(&d_compacted_size, sizeof(uint64_t));

	OptixAccelEmitDesc _emit_desc;
	_emit_desc.type	= OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
	_emit_desc.result	= (CUdeviceptr) d_compacted_size;

	void * d_temp_buffer;
	cudaMalloc(&d_temp_buffer, _gas_buffer_size.tempSizeInBytes);

	void * d_output_buffer;
	cudaMalloc(&d_output_buffer, _gas_buffer_size.outputSizeInBytes);


	optixAccelBuild(	_context,
						0,	// stream
						&_accel_opt,
						_meshes.data(),
						size(_meshes),
						(CUdeviceptr) d_temp_buffer,
						_gas_buffer_size.tempSizeInBytes,
						(CUdeviceptr) d_output_buffer,
						_gas_buffer_size.outputSizeInBytes,
						&_as_handle,
						&_emit_desc,
						1
						);

	cudaDeviceSynchronize();

	uint64_t compacted_size;
	cudaMemcpy(&compacted_size, d_compacted_size, sizeof(uint64_t), cudaMemcpyDeviceToHost);

	cudaMalloc(&as_buffer, compacted_size);

	optixAccelCompact(	_context,
						0,	// stream
						_as_handle,
						(CUdeviceptr) as_buffer,
						compacted_size,
						&_as_handle
						);

	cudaDeviceSynchronize();

	cudaFree(d_output_buffer);
	cudaFree(d_temp_buffer);
	cudaFree(d_compacted_size);

	return _as_handle;
}

void optix::add_mesh(OptixBuildInput & _mesh, CUdeviceptr & d_vertex_ptr, uint32_t & _trig_flags, const che * mesh, const mat4 & model_mat)
{
	che * d_m = new che_cuda(mesh);
	d_mesh.push_back(d_m);

	float * d_model_mat = nullptr;
	cudaMalloc(&d_model_mat, sizeof(model_mat));
	cudaMemcpy(d_model_mat, &model_mat, sizeof(model_mat), cudaMemcpyHostToDevice);

	d_vertex_ptr = (CUdeviceptr) &d_m->point(0);

	_mesh = {};
	_mesh.type = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;

	_mesh.triangleArray.vertexFormat		= OPTIX_VERTEX_FORMAT_FLOAT3;
	_mesh.triangleArray.vertexStrideInBytes	= 3 * sizeof(float);
	_mesh.triangleArray.numVertices			= d_m->n_vertices;
	_mesh.triangleArray.vertexBuffers		= &d_vertex_ptr;

	_mesh.triangleArray.indexFormat			= OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
	_mesh.triangleArray.indexStrideInBytes	= 3 * sizeof(index_t);
	_mesh.triangleArray.numIndexTriplets	= d_m->n_trigs;
	_mesh.triangleArray.indexBuffer			= (CUdeviceptr) d_m->trigs_ptr();

	_mesh.triangleArray.transformFormat		= OPTIX_TRANSFORM_FORMAT_MATRIX_FLOAT12;
	_mesh.triangleArray.preTransform		= (CUdeviceptr) d_model_mat;

	_trig_flags = 0;

	_mesh.triangleArray.flags						= &_trig_flags;
	_mesh.triangleArray.numSbtRecords				= 1;
	_mesh.triangleArray.sbtIndexOffsetBuffer		= 0;
	_mesh.triangleArray.sbtIndexOffsetSizeInBytes	= 0;
	_mesh.triangleArray.sbtIndexOffsetStrideInBytes	= 0;

	if(mesh->is_scene())
	{
		scene * sc = (scene *) mesh;
		cudaMalloc(&params.sc.materials, size(sc->materials) * sizeof(scene::material));
		cudaMalloc(&params.sc.textures, size(sc->textures) * sizeof(scene::texture));
		cudaMalloc(&params.sc.trig_mat, mesh->n_vertices / 3 * sizeof(index_t));
		cudaMalloc(&params.sc.texcoords, mesh->n_vertices * sizeof(vec2));

		std::vector<scene::texture> textures = sc->textures;
		for(scene::texture & tex: textures)
		{
			unsigned char * h_data = tex.data;
			cudaMalloc(&tex.data, tex.width * tex.height * tex.spectrum);
			cudaMemcpy(tex.data, h_data, tex.width * tex.height * tex.spectrum, cudaMemcpyHostToDevice);
			tex_data.push_back(tex.data);
		}

		gproshan_error_var(size(textures));
		cudaMemcpy(params.sc.materials, sc->materials.data(), size(sc->materials) * sizeof(scene::material), cudaMemcpyHostToDevice);
		cudaMemcpy(params.sc.textures, textures.data(), size(textures) * sizeof(scene::texture), cudaMemcpyHostToDevice);
		cudaMemcpy(params.sc.trig_mat, sc->trig_mat, mesh->n_vertices / 3 * sizeof(index_t), cudaMemcpyHostToDevice);
		cudaMemcpy(params.sc.texcoords, sc->texcoords, mesh->n_vertices * sizeof(vec2), cudaMemcpyHostToDevice);
	}
}


} // namespace gproshan

#endif // GPROSHAN_OPTIX

