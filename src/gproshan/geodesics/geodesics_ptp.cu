#include <gproshan/geodesics/geodesics_ptp.h>

#include <gproshan/mesh/che_cuda.h>

#include <cstdio>
#include <fstream>
#include <cassert>
#include <cublas_v2.h>


// geometry processing and shape analysis framework
namespace gproshan {


double parallel_toplesets_propagation_gpu(	const ptp_out_t & ptp_out,
											const che * mesh,
											const std::vector<index_t> & sources,
											const toplesets & tps,
											const bool coalescence,
											const bool set_inf,
											const f_ptp<float> & fun
											)
{
	cudaDeviceReset();

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);


	const coalescence_ptp inv(coalescence ? mesh : nullptr, tps);
	const size_t n_vertices = coalescence ? inv->n_vertices : mesh->n_vertices;

	const che_cuda d_mesh(coalescence ? inv : mesh, {false, false, false});

	gproshan_error_var(coalescence);
	gproshan_error_var(n_vertices == mesh->n_vertices);
	gproshan_error_var(n_vertices == tps.size());

	float * h_dist = new float[n_vertices];
	index_t * h_clusters = ptp_out.clusters ? new index_t[n_vertices] : nullptr;

	float * d_dist[3] = {};
	index_t * d_clusters[3] = {};
	index_t * d_sorted = nullptr;
	index_t * d_inv = nullptr;

	cudaMalloc(&d_dist[0], sizeof(float) * n_vertices);
	cudaMalloc(&d_dist[1], sizeof(float) * n_vertices);
	d_dist[2] = h_dist;

	if(h_clusters)
	{
		cudaMalloc(&d_clusters[0], sizeof(index_t) * n_vertices);
		cudaMalloc(&d_clusters[1], sizeof(index_t) * n_vertices);
		d_clusters[2] = h_clusters;
	}

	if(!coalescence)
	{
		cudaMalloc(&d_sorted, sizeof(index_t) * n_vertices);
		cudaMemcpy(d_sorted, tps.sorted, sizeof(index_t) * std::size(tps), cudaMemcpyHostToDevice);

		cudaMalloc(&d_inv, sizeof(index_t) * mesh->n_vertices);
		cudaMemcpy(d_inv, tps.sorted, sizeof(index_t) * mesh->n_vertices, cudaMemcpyHostToDevice);
	}

	if(set_inf)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
			h_dist[v] = INFINITY;
	}

	const index_t i = run_ptp(d_mesh, sources, tps.splits, d_dist, d_clusters, d_sorted, d_inv, fun);

	cudaMemcpy(h_dist, d_dist[i], sizeof(float) * n_vertices, cudaMemcpyDeviceToHost);

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		ptp_out.dist[tps.sorted[v]] = h_dist[v];

	delete [] h_dist;

	if(h_clusters)
	{
		cudaMemcpy(h_clusters, d_clusters[i], sizeof(index_t) * n_vertices, cudaMemcpyDeviceToHost);

		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
			ptp_out.clusters[tps.sorted[v]] = h_clusters[v];

		delete [] h_clusters;
	}

	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cudaFree(d_clusters[0]);
	cudaFree(d_clusters[1]);
	cudaFree(d_sorted);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	float time;
	cudaEventElapsedTime(&time, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return time / 1000;
}

double farthest_point_sampling_ptp_gpu(che * mesh, std::vector<index_t> & samples, size_t n, float radio)
{
	const size_t n_vertices = mesh->n_vertices;

	cudaDeviceReset();

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	const che_cuda d_mesh(mesh, {false, false, false});

	float * h_dist = new float[n_vertices];

	float * d_dist[3] = {};
	index_t * d_clusters[3] = {};
	index_t * d_sorted = nullptr;
	index_t * d_inv = nullptr;

	cudaMalloc(&d_dist[0], sizeof(float) * n_vertices);
	cudaMalloc(&d_dist[1], sizeof(float) * n_vertices);
	cudaMalloc(&d_sorted, sizeof(index_t) * n_vertices);
	cudaMalloc(&d_inv, sizeof(index_t) * mesh->n_vertices);
	d_dist[2] = h_dist;

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		h_dist[v] = INFINITY;

	toplesets tps(mesh, samples);

	cublasHandle_t handle;
	cublasCreate(&handle);

	if(n >= n_vertices) n = n_vertices >> 2;

	n -= size(samples);
	samples.reserve(n);

	int farthest;
	float max_dist = INFINITY;
	while(n-- && radio < max_dist)
	{
		cudaMemcpy(d_sorted, tps.sorted, sizeof(index_t) * tps.size(), cudaMemcpyHostToDevice);
		const index_t i = run_ptp(d_mesh, samples, tps.splits, d_dist, d_clusters, d_sorted, d_inv);

		// 1 indexing
		cublasIsamax(handle, mesh->n_vertices, d_dist[i], 1, &farthest);

		if(radio > 0 || !n)
			cudaMemcpy(&max_dist, d_dist[i] + farthest - 1, sizeof(float), cudaMemcpyDeviceToHost);

		samples.push_back(tps.sorted[farthest - 1]);
		tps.reset(mesh, samples);
	}

	cublasDestroy(handle);

	delete [] h_dist;

	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cudaFree(d_sorted);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	float time;
	cudaEventElapsedTime(&time, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return time / 1000;
}

__global__
void relax_ptp(const che * mesh, float * new_dist, float * old_dist, index_t * new_clusters, index_t * old_clusters, const index_t start, const index_t end, const index_t * sorted, const index_t * inv)
{
	index_t i = blockDim.x * blockIdx.x + threadIdx.x + start;
	if(i >= end) return;

	relax_ptp(mesh, sorted, inv, i, new_dist, old_dist, new_clusters, old_clusters);
}

__global__
void relative_error(unsigned int * g_count, const float * new_dist, const float * old_dist, const index_t start, const index_t end)
{
	const index_t i = blockDim.x * blockIdx.x + threadIdx.x + start;
	if(i >= end) return;

	__shared__ unsigned int count;
	if(!threadIdx.x)
		count = 0;

	if(!i) *g_count = 0;

	__syncthreads();

	atomicInc(&count, fabsf(new_dist[i] - old_dist[i]) / old_dist[i] < PTP_TOL);

	__syncthreads();

	if(!threadIdx.x)
		atomicInc(g_count, count);
}


} // namespace gproshan

