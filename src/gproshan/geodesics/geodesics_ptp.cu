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
											const toplesets_t & toplesets,
											const bool coalescence,
											const bool set_inf,
											const f_ptp<float> & fun
											)
{
	const size_t n_vertices = mesh->n_vertices;

	che * h_mesh = nullptr;
	index_t * inv = nullptr;
	if(coalescence)
	{
		h_mesh = new che(*mesh, toplesets.index, {false, false, false});
		inv = new index_t[n_vertices];

		#pragma omp parallel for
		for(index_t i = 0; i < toplesets.limits.back(); ++i)
			inv[toplesets.index[i]] = i;
	}


	cudaDeviceReset();

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	che_cuda d_mesh(h_mesh ? h_mesh : mesh, {false, false, false});

	float * h_dist = coalescence ? new float[n_vertices] : ptp_out.dist;
	index_t * h_clusters = coalescence && ptp_out.clusters ? new index_t[n_vertices]
															: ptp_out.clusters;

	float * d_error = nullptr;
	float * d_dist[3] = {};
	index_t * d_clusters[3] = {};
	index_t * d_sorted = nullptr;

	cudaMalloc(&d_error, sizeof(float) * n_vertices);
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
	}

	if(set_inf)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
			h_dist[v] = INFINITY;
	}

	const index_t i = run_ptp(	d_mesh, sources, toplesets.limits, d_error, d_dist, d_clusters,
								coalescence ? inv : toplesets.index, d_sorted,
								fun);

	cudaMemcpy(h_dist, d_dist[i], sizeof(float) * n_vertices, cudaMemcpyDeviceToHost);

	if(coalescence)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
			ptp_out.dist[v] = h_dist[inv[v]];

		delete [] h_dist;
	}

	if(h_clusters)
	{
		cudaMemcpy(h_clusters, d_clusters[i], sizeof(index_t) * n_vertices, cudaMemcpyDeviceToHost);

		if(coalescence)
		{
			#pragma omp parallel for
			for(index_t v = 0; v < n_vertices; ++v)
				ptp_out.clusters[v] = h_clusters[inv[v]];

			delete [] h_clusters;
		}
	}

	cudaFree(d_error);
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cudaFree(d_clusters[0]);
	cudaFree(d_clusters[1]);
	cudaFree(d_sorted);

	delete [] inv;
	delete h_mesh;

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	float time;
	cudaEventElapsedTime(&time, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return time / 1000;
}

float farthest_point_sampling_ptp_gpu(che * mesh, std::vector<index_t> & samples, double & time_fps, size_t n, float radio)
{
	const size_t n_vertices = mesh->n_vertices;

	cudaDeviceReset();

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	che_cuda d_mesh(mesh, {false, false, false});

	float * h_dist = new float[n_vertices];

	float * d_error = nullptr;
	float * d_dist[3] = {};
	index_t * d_clusters[3] = {};
	index_t * d_sorted = nullptr;

	cudaMalloc(&d_error, sizeof(float) * n_vertices);
	cudaMalloc(&d_dist[0], sizeof(float) * n_vertices);
	cudaMalloc(&d_dist[1], sizeof(float) * n_vertices);
	cudaMalloc(&d_sorted, sizeof(index_t) * n_vertices);
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
		const index_t i = run_ptp(d_mesh, samples, tps.splits, d_error, d_dist, d_clusters, tps.sorted, d_sorted);

		// 1 indexing
		cublasIsamax(handle, mesh->n_vertices, d_dist[i], 1, &farthest);

		if(radio > 0 || !n)
			cudaMemcpy(&max_dist, d_dist[i] + farthest - 1, sizeof(float), cudaMemcpyDeviceToHost);

		samples.push_back(farthest - 1);
		tps.reset(mesh, samples);
	}

	cublasDestroy(handle);

	delete [] h_dist;

	cudaFree(d_error);
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cudaFree(d_sorted);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	float time;
	cudaEventElapsedTime(&time, start, stop);
	time_fps = time / 1000;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return max_dist;
}

__global__
void relax_ptp(const che * mesh, float * new_dist, float * old_dist, index_t * new_clusters, index_t * old_clusters, const index_t start, const index_t end, const index_t * sorted)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;

	if(v < end)
		relax_ptp(mesh, new_dist, old_dist, new_clusters, old_clusters, sorted ? sorted[v] : v);
}

__global__
void relative_error(float * error, const float * new_dist, const float * old_dist, const index_t start, const index_t end, const index_t * sorted)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;

	if(v < end)
	{
		v = sorted ? sorted[v] : v;
		error[v] = fabsf(new_dist[v] - old_dist[v]) / old_dist[v];
	}
}

__host_device__
bool is_ok::operator()(const float val) const
{
	return val < PTP_TOL;
}

__host_device__
bool is_ok::operator()(const index_t i) const
{
	return error[i] < PTP_TOL;
}


} // namespace gproshan

