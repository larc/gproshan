#include "mesh/che.cuh"


// geometry processing and shape analysis framework
namespace gproshan {


__host__ __device__
index_t cu_trig(index_t he)
{
	if(he == NIL) return NIL;
	return he / che::mtrig;
}

__host__ __device__
index_t cu_next(index_t he)
{
	if(he == NIL) return NIL;
	return che::mtrig * cu_trig(he) + (he + 1) % che::mtrig;
}

__host__ __device__
index_t cu_prev(index_t he)
{
	if(he == NIL) return NIL;
	return che::mtrig * cu_trig(he) + (he + che::mtrig - 1) % che::mtrig;
}

void cuda_create_CHE(CHE * h_che, CHE *& dd_che, CHE *& d_che, const bool & normal)
{
	dd_che = (CHE *) malloc(sizeof(CHE));
	memcpy(dd_che, h_che, sizeof(CHE));

	cudaMalloc(&dd_che->GT, sizeof(vertex_cu) * h_che->n_vertices);
	cudaMemcpy(dd_che->GT, h_che->GT, sizeof(vertex_cu) * h_che->n_vertices, cudaMemcpyHostToDevice);
	
	if(normal)
	{
		cudaMalloc(&dd_che->VN, sizeof(vertex_cu) * h_che->n_vertices);
		cudaMemcpy(dd_che->VN, h_che->VN, sizeof(vertex_cu) * h_che->n_vertices, cudaMemcpyHostToDevice);
	}

	cudaMalloc(&dd_che->VT, sizeof(index_t) * h_che->n_half_edges);
	cudaMemcpy(dd_che->VT, h_che->VT, sizeof(index_t) * h_che->n_half_edges, cudaMemcpyHostToDevice);

	cudaMalloc(&dd_che->OT, sizeof(index_t) * h_che->n_half_edges);
	cudaMemcpy(dd_che->OT, h_che->OT, sizeof(index_t) * h_che->n_half_edges, cudaMemcpyHostToDevice);

	cudaMalloc(&dd_che->EVT, sizeof(index_t) * h_che->n_vertices);
	cudaMemcpy(dd_che->EVT, h_che->EVT, sizeof(index_t) * h_che->n_vertices, cudaMemcpyHostToDevice);

	cudaMalloc(&d_che, sizeof(CHE));
	cudaMemcpy(d_che, dd_che, sizeof(CHE), cudaMemcpyHostToDevice);
}

void cuda_free_CHE(CHE *& dd_che, CHE *& d_che)
{
	if(dd_che->GT) cudaFree(dd_che->GT);
	if(dd_che->VN) cudaFree(dd_che->VN);
	if(dd_che->VT) cudaFree(dd_che->VT);
	if(dd_che->OT) cudaFree(dd_che->OT);
	if(dd_che->EVT) cudaFree(dd_che->EVT);

	free(dd_che);
	cudaFree(d_che);
}


} // namespace gproshan

