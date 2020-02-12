#include "rt_optix.h"

#ifdef GPROSHAN_OPTIX

#include <iostream>
#include <random>
#include <cstring>

#include <optix_function_table_definition.h>

// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


optix::optix(const std::vector<che *> & meshes)
{
	optixInit();
}

optix::~optix()
{
}

OptixTraversableHandle optix::build_as(const std::vector<che *> & meshes)
{
	std::vector<OptixBuildInput> optix_meshes(meshes.size());
	std::vector<uint32_t> optix_trig_flags(meshes.size());

	for(index_t i = 0; i < meshes.size(); i++)
		add_mesh(optix_meshes[i], optix_trig_flags[i], meshes[i]);

	OptixTraversableHandle optix_as_handle = {};
	
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
	optix_mesh.triangleArray.numVertices			= (int) mesh->n_vertices(); 
	optix_mesh.triangleArray.vertexBuffers			= (CUdeviceptr *) &d_vertex;

	optix_mesh.triangleArray.indexFormat			= OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
	optix_mesh.triangleArray.indexStrideInBytes		= 3 * sizeof(index_t);
	optix_mesh.triangleArray.numIndexTriplets		= (int) mesh->n_faces();
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

