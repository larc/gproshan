#include "rt_optix.h"

#ifdef GPROSHAN_OPTIX

#include <iostream>
#include <random>
#include <cstring>

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
	
	cudaStreamCreate(&stream);

	cudaGetDeviceProperties(&device_prop, 0);	// device id = 0

	cuCtxGetCurrent(&cuda_context);

	optixDeviceContextCreate(cuda_context, 0, &optix_context);
	optixDeviceContextSetLogCallback(optix_context, optix_log, nullptr, 4);

	build_as(meshes);
}

optix::~optix()
{
}

OptixTraversableHandle optix::build_as(const std::vector<che *> & meshes)
{
	OptixTraversableHandle optix_as_handle = {};
	
	std::vector<OptixBuildInput> optix_meshes(meshes.size());
	std::vector<uint32_t> optix_trig_flags(meshes.size());

	for(index_t i = 0; i < meshes.size(); i++)
		add_mesh(optix_meshes[i], optix_trig_flags[i], meshes[i]);

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

void optix::add_mesh(OptixBuildInput & optix_mesh, uint32_t & optix_trig_flags, const che * mesh)
{
	void * d_vertex;
	void * d_index;
	
	cudaMalloc(&d_vertex, mesh->n_vertices() * sizeof(vertex));
	cudaMemcpy(d_vertex, &mesh->gt(0), mesh->n_vertices() * sizeof(vertex), cudaMemcpyHostToDevice);
	
	cudaMalloc(&d_index, mesh->n_half_edges() * sizeof(index_t));
	cudaMemcpy(d_index, &mesh->vt(0), mesh->n_half_edges() * sizeof(index_t), cudaMemcpyHostToDevice);

	optix_mesh = {};
	optix_mesh.type = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;

#ifdef SINGLE_P
	optix_mesh.triangleArray.vertexFormat			= OPTIX_VERTEX_FORMAT_FLOAT3;
#else
	optix_mesh.triangleArray.vertexFormat			= OPTIX_VERTEX_FORMAT_DOUBLE3; //ERROR! 
#endif // SINGLE_P
	optix_mesh.triangleArray.vertexStrideInBytes	= sizeof(vertex); 
	optix_mesh.triangleArray.numVertices			= mesh->n_vertices(); 
	optix_mesh.triangleArray.vertexBuffers			= (CUdeviceptr *) &d_vertex;

	optix_mesh.triangleArray.indexFormat			= OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
	optix_mesh.triangleArray.indexStrideInBytes		= 3 * sizeof(index_t);
	optix_mesh.triangleArray.numIndexTriplets		= mesh->n_faces();
	optix_mesh.triangleArray.indexBuffer			= (CUdeviceptr) d_index;

	optix_trig_flags = 0 ;
	
	optix_mesh.triangleArray.flags							= &optix_trig_flags;
	optix_mesh.triangleArray.numSbtRecords					= 1;
	optix_mesh.triangleArray.sbtIndexOffsetBuffer			= 0; 
	optix_mesh.triangleArray.sbtIndexOffsetSizeInBytes		= 0; 
	optix_mesh.triangleArray.sbtIndexOffsetStrideInBytes	= 0; 
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

