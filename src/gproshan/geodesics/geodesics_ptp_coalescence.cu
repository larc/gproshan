#include <gproshan/geodesics/geodesics_ptp_coalescence.cuh>

#include <gproshan/mesh/che_off.h>

#include <cstdio>
#include <fstream>
#include <cassert>
#include <cublas_v2.h>

#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>


// geometry processing and shape analysis framework
namespace gproshan {


double parallel_toplesets_propagation_coalescence_gpu(const ptp_out_t & ptp_out, const che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets, const bool & set_inf)
{
	index_t * inv = nullptr;
	che * coalescence_mesh = ptp_coalescence(inv, mesh, toplesets);

	// ------------------------------------------------------

	cudaDeviceReset();

	float time;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	// BEGIN PTP

	CHE * h_mesh = new CHE(coalescence_mesh);
	CHE * dd_mesh, * d_mesh;
	cuda_create_CHE(h_mesh, dd_mesh, d_mesh);

	real_t * h_dist = new real_t[h_mesh->n_vertices];

	if(set_inf)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < h_mesh->n_vertices; ++v)
			h_dist[v] = INFINITY;
	}

	real_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(real_t) * h_mesh->n_vertices);
	cudaMalloc(&d_dist[1], sizeof(real_t) * h_mesh->n_vertices);

	real_t * d_error;
	cudaMalloc(&d_error, sizeof(real_t) * h_mesh->n_vertices);

	index_t d;
	if(ptp_out.clusters)
	{
		index_t * h_clusters = new index_t[h_mesh->n_vertices];
		index_t * d_clusters[2] = {nullptr, nullptr};

		cudaMalloc(&d_clusters[0], sizeof(index_t) * h_mesh->n_vertices);
		cudaMalloc(&d_clusters[1], sizeof(index_t) * h_mesh->n_vertices);

		d = run_ptp_coalescence_gpu(d_mesh, h_mesh->n_vertices, h_dist, d_dist, sources, {toplesets.limits, inv}, d_error, h_clusters, d_clusters);
		cudaMemcpy(h_clusters, d_clusters[d], sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);

		#pragma omp parallel for
		for(index_t i = 0; i < h_mesh->n_vertices; ++i)
			ptp_out.clusters[toplesets.index[i]] = h_clusters[i];

		cudaFree(d_clusters[0]);
		cudaFree(d_clusters[1]);

		delete [] h_clusters;
	}
	else d = run_ptp_coalescence_gpu(d_mesh, h_mesh->n_vertices, h_dist, d_dist, sources, {toplesets.limits, inv}, d_error);

	cudaMemcpy(h_dist, d_dist[d], sizeof(real_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);

	cudaFree(d_error);
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cuda_free_CHE(dd_mesh, d_mesh);

	delete coalescence_mesh;
	delete [] inv;

	#pragma omp parallel for
	for(index_t i = 0; i < toplesets.limits.back(); ++i)
		ptp_out.dist[toplesets.index[i]] = h_dist[i];

	delete [] h_dist;

	// END PTP

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return time / 1000;
}

index_t run_ptp_coalescence_gpu(CHE * d_mesh, const index_t & n_vertices, real_t * h_dist, real_t ** d_dist, const std::vector<index_t> & sources, const toplesets_t & inv, real_t * d_error, index_t * h_clusters, index_t ** d_clusters)
{
	for(index_t i = 0; i < sources.size(); ++i)
		h_dist[inv.index[sources[i]]] = 0;

	cudaMemcpy(d_dist[0], h_dist, sizeof(real_t) * n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dist[1], h_dist, sizeof(real_t) * n_vertices, cudaMemcpyHostToDevice);

	if(h_clusters)
	{
		assert(d_clusters[0]);

		for(index_t i = 0; i < sources.size(); ++i)
			h_clusters[inv.index[sources[i]]] = i + 1;

		cudaMemcpy(d_clusters[0], h_clusters, sizeof(index_t) * n_vertices, cudaMemcpyHostToDevice);
		cudaMemcpy(d_clusters[1], h_clusters, sizeof(index_t) * n_vertices, cudaMemcpyHostToDevice);
	}

	index_t d = 0;
	index_t start, end, n_cond;
	index_t i = 1, j = 2;

	// maximum number of iterations
	index_t iter = 0;
	index_t max_iter = inv.limits.size() << 1;

	while(i < j && iter++ < max_iter)
	{
		if(i < (j >> 1)) i = (j >> 1); // K/2 limit band size

		start = inv.limits[i];
		end = inv.limits[j];
		n_cond = inv.limits[i + 1] - start;

		h_clusters ? relax_ptp <<< NB(end - start), NT >>> (d_mesh, d_dist[!d], d_dist[d], d_clusters[!d], d_clusters[d], start, end)
					: relax_ptp <<< NB(end - start), NT >>> (d_mesh, d_dist[!d], d_dist[d], nullptr, nullptr, start, end);

		cudaDeviceSynchronize();

		relative_error <<< NB(n_cond), NT >>>(d_error, d_dist[!d], d_dist[d], start, start + n_cond);
		cudaDeviceSynchronize();

		if(n_cond == thrust::count_if(thrust::device, d_error + start, d_error + start + n_cond, is_ok()))
			++i;

		if(j < inv.limits.size() - 1) ++j;

		d = !d;
	}

	return d;
}


} // namespace gproshan

