#include <gproshan/mesh/che.cuh>


// geometry processing and shape analysis framework
namespace gproshan {


void cuda_create_CHE(const che * h_che, che *& dd_che, che *& d_che, const bool normal, const bool color)
{
	dd_che = new che;
	dd_che->n_vertices = h_che->n_vertices;
	dd_che->n_trigs = h_che->n_trigs;
	dd_che->n_half_edges = h_che->n_half_edges;

	cudaMalloc(&dd_che->GT, sizeof(vertex) * h_che->n_vertices);
	cudaMemcpy(dd_che->GT, h_che->GT, sizeof(vertex) * h_che->n_vertices, cudaMemcpyHostToDevice);

	if(normal)
	{
		cudaMalloc(&dd_che->VN, sizeof(vertex) * h_che->n_vertices);
		cudaMemcpy(dd_che->VN, h_che->VN, sizeof(vertex) * h_che->n_vertices, cudaMemcpyHostToDevice);
	}

	if(color)
	{
		cudaMalloc(&dd_che->VC, sizeof(che::rgb_t) * h_che->n_vertices);
		cudaMemcpy(dd_che->VC, h_che->VC, sizeof(che::rgb_t) * h_che->n_vertices, cudaMemcpyHostToDevice);

		cudaMalloc(&dd_che->VHC, sizeof(real_t) * h_che->n_vertices);
		cudaMemcpy(dd_che->VHC, h_che->VHC, sizeof(real_t) * h_che->n_vertices, cudaMemcpyHostToDevice);
	}

	cudaMalloc(&dd_che->VT, sizeof(index_t) * h_che->n_half_edges);
	cudaMemcpy(dd_che->VT, h_che->VT, sizeof(index_t) * h_che->n_half_edges, cudaMemcpyHostToDevice);

	cudaMalloc(&dd_che->OT, sizeof(index_t) * h_che->n_half_edges);
	cudaMemcpy(dd_che->OT, h_che->OT, sizeof(index_t) * h_che->n_half_edges, cudaMemcpyHostToDevice);

	cudaMalloc(&dd_che->EVT, sizeof(index_t) * h_che->n_vertices);
	cudaMemcpy(dd_che->EVT, h_che->EVT, sizeof(index_t) * h_che->n_vertices, cudaMemcpyHostToDevice);

	cudaMalloc(&d_che, sizeof(che));
	cudaMemcpy(d_che, dd_che, sizeof(che), cudaMemcpyHostToDevice);
}

void cuda_free_CHE(che *& dd_che, che *& d_che)
{
	cudaFree(dd_che->GT);
	cudaFree(dd_che->EVT);
	cudaFree(dd_che->VT);
	cudaFree(dd_che->OT);
	cudaFree(dd_che->VN);
	cudaFree(dd_che->VC);
	cudaFree(dd_che->VHC);

	//delete dd_che;
	cudaFree(d_che);
}


} // namespace gproshan

