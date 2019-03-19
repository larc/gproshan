#include "test_geodesics_ptp.cuh"
#include "test_geodesics_ptp.h"

#include "geodesics_ptp.cuh"
#include "geodesics_ptp.h"

#include <fstream>
#include <cublas_v2.h>

distance_t * iter_error_parallel_toplesets_propagation_gpu(che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, const distance_t * exact_dist, double & time_ptp)
{
	cudaDeviceReset();

	float time;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	// BEGIN PTP

	CHE * h_mesh = new CHE(mesh);
	CHE * dd_mesh, * d_mesh;
	cuda_create_CHE(h_mesh, dd_mesh, d_mesh);

	distance_t * h_dist = new distance_t[h_mesh->n_vertices];

	distance_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(distance_t) * h_mesh->n_vertices);
	cudaMalloc(&d_dist[1], sizeof(distance_t) * h_mesh->n_vertices);

	index_t * d_sorted;
	cudaMalloc(&d_sorted, sizeof(index_t) * h_mesh->n_vertices);

	distance_t * error = iter_error_run_ptp_gpu(d_mesh, h_mesh->n_vertices, h_dist, d_dist, sources, limits, sorted_index, d_sorted, exact_dist);
	
	delete [] h_dist;
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cudaFree(d_sorted);
	cuda_free_CHE(dd_mesh, d_mesh);

	// END PTP

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	time_ptp = time / 1000;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	
	return error;
}

/// Return an array of time in seconds.
double * times_farthest_point_sampling_ptp_gpu(che * mesh, vector<index_t> & samples, size_t n, distance_t radio)
{
	cudaDeviceReset();

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// BEGIN FPS PTP

	CHE * h_mesh = new CHE(mesh);
	CHE * dd_mesh, * d_mesh;
	cuda_create_CHE(h_mesh, dd_mesh, d_mesh);

	distance_t * h_dist = new distance_t[h_mesh->n_vertices];

	distance_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(distance_t) * h_mesh->n_vertices);
	cudaMalloc(&d_dist[1], sizeof(distance_t) * h_mesh->n_vertices);

	distance_t * d_error;
	cudaMalloc(&d_error, sizeof(distance_t) * h_mesh->n_vertices);

	index_t * d_sorted;
	cudaMalloc(&d_sorted, sizeof(index_t) * h_mesh->n_vertices);

	vector<index_t> limits;
	index_t * toplesets = new index_t[h_mesh->n_vertices];
	index_t * sorted_index = new index_t[h_mesh->n_vertices];

	cublasHandle_t handle;
	cublasCreate(&handle);

	if(n >= h_mesh->n_vertices) n = h_mesh->n_vertices >> 1;

	double * times = new double[n + 1];

	n -= samples.size();
	samples.reserve(n);

	float time_fps;
	index_t d;
	int f;
	distance_t max_dist = INFINITY;
	while(n-- && max_dist > radio)
	{
		cudaEventRecord(start, 0);
		
		limits.clear();
		mesh->compute_toplesets(toplesets, sorted_index, limits, samples);

		d = run_ptp_gpu(d_mesh, h_mesh->n_vertices, h_dist, d_dist, samples, limits, sorted_index, d_sorted, d_error);

		// 1 indexing
		#ifdef SINGLE_P
			cublasIsamax(handle, mesh->n_vertices(), d_dist[d], 1, &f);
		#else
			cublasIdamax(handle, mesh->n_vertices(), d_dist[d], 1, &f);
		#endif
		
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time_fps, start, stop);

		times[samples.size()] = time_fps / 1000;

		if(radio > 0 || !n)
			cudaMemcpy(&max_dist, d_dist[d] + f - 1, sizeof(distance_t), cudaMemcpyDeviceToHost);

		samples.push_back(f - 1);
	}

	cublasDestroy(handle);

	delete [] h_dist;
	delete [] toplesets;
	delete [] sorted_index;
	
	cudaFree(d_error);
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cudaFree(d_sorted);
	cuda_free_CHE(dd_mesh, d_mesh);

	// END FPS PTP

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return times;
}

distance_t * iter_error_run_ptp_gpu(CHE * d_mesh, const index_t & n_vertices, distance_t * h_dist, distance_t ** d_dist, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * h_sorted, index_t * d_sorted, const distance_t * exact_dist)
{
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; v++)
		h_dist[v] = INFINITY;

	for(index_t i = 0; i < sources.size(); i++)
		h_dist[sources[i]] = 0;

	cudaMemcpy(d_dist[0], h_dist, sizeof(distance_t) * n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dist[1], h_dist, sizeof(distance_t) * n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_sorted, h_sorted, sizeof(index_t) * n_vertices, cudaMemcpyHostToDevice);

	index_t d = 0, e = 0;
	index_t start, end;
	index_t iter = iterations(limits);

	distance_t * dist_error = new distance_t[iter - limits.size()];	

	for(index_t i = 2; i < iter; i++)
	{
		start = start_v(i, limits);
		end = end_v(i, limits);

		if(end == start) break;

		relax_ptp <<< NB(end - start), NT >>> (d_mesh, d_dist[!d], d_dist[d], d_sorted, end, start);
		cudaMemcpy(h_dist, d_dist[!d], sizeof(distance_t) * n_vertices, cudaMemcpyDeviceToHost);
		
		// calculating iteration error
		if(i >= limits.size())
			dist_error[e++] = compute_error(h_dist, exact_dist, n_vertices, sources.size());

		d = !d;
	}

	return dist_error;
}

