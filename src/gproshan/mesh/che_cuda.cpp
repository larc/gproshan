#include <gproshan/mesh/che_cuda.h>

#include <cstring>
#include <cstdio>
#include <cassert>

#include <cuda_runtime.h>


// geometry processing and shape analysis framework
namespace gproshan {


che_cuda::che_cuda(const che * mesh, const che::options & opts)
{
	if(!mesh) return;
	if(!mesh->n_vertices) return;

	n_vertices = mesh->n_vertices;
	n_trigs = mesh->n_trigs;
	n_half_edges = mesh->n_half_edges;

	cudaMalloc(&GT, sizeof(vertex) * n_vertices);
	cudaMemcpy(GT, mesh->GT, sizeof(vertex) * n_vertices, cudaMemcpyHostToDevice);

	cudaMalloc(&EVT, sizeof(index_t) * mesh->n_vertices);
	cudaMemcpy(EVT, mesh->EVT, sizeof(index_t) * mesh->n_vertices, cudaMemcpyHostToDevice);

	cudaMalloc(&VT, sizeof(index_t) * n_half_edges);
	cudaMemcpy(VT, mesh->VT, sizeof(index_t) * n_half_edges, cudaMemcpyHostToDevice);

	cudaMalloc(&OT, sizeof(index_t) * n_half_edges);
	cudaMemcpy(OT, mesh->OT, sizeof(index_t) * n_half_edges, cudaMemcpyHostToDevice);

	if(opts.normals)
	{
		cudaMalloc(&VN, sizeof(vertex) * n_vertices);
		cudaMemcpy(VN, mesh->VN, sizeof(vertex) * n_vertices, cudaMemcpyHostToDevice);
	}

	if(opts.colors)
	{
		cudaMalloc(&VC, sizeof(che::rgb_t) * n_vertices);
		cudaMemcpy(VC, mesh->VC, sizeof(che::rgb_t) * n_vertices, cudaMemcpyHostToDevice);

		cudaMalloc(&VHC, sizeof(real_t) * n_vertices);
		cudaMemcpy(VHC, mesh->VHC, sizeof(real_t) * n_vertices, cudaMemcpyHostToDevice);
	}

	cudaMalloc(&device, sizeof(che));
	cudaMemcpy(device, this, sizeof(che), cudaMemcpyHostToDevice);
}

che_cuda::~che_cuda()
{
	gproshan_error_var("che_cuda::free: " + filename);

	cudaFree(GT);	GT	= nullptr;
	cudaFree(EVT);	EVT	= nullptr;
	cudaFree(VT);	VT	= nullptr;
	cudaFree(OT);	OT	= nullptr;
	cudaFree(VN);	VN	= nullptr;
	cudaFree(VC);	VC	= nullptr;
	cudaFree(VHC);	VHC	= nullptr;

	cudaFree(device);
}

che_cuda::operator const che * () const
{
	return device;
}


} // namespace gproshan

