#include "test_geodesics_ptp_coalescence.cuh"

#include "geodesics_ptp_coalescence.cuh"
#include "geodesics_ptp.h"

#include "che_off.h"

#include <fstream>
#include <omp.h>
#include <cublas_v2.h>

distance_t * iter_error_parallel_toplesets_propagation_coalescence_gpu(che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, const distance_t * exact_dist, double & time_ptp)
{
	// sort data by levels, must be improve the coalescence

	vertex * V = new vertex[mesh->n_vertices()];
	index_t * F = new index_t[mesh->n_faces() * che::P];
	index_t * inv = new index_t[mesh->n_vertices()];
	
	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices(); i++)
	{
		V[i] = mesh->gt(sorted_index[i]);
		inv[sorted_index[i]] = i;
	}

	#pragma omp parallel for
	for(index_t he = 0; he < mesh->n_half_edges(); he++)
		F[he] = inv[mesh->vt(he)];

	mesh = new che_off(V, mesh->n_vertices(), F, mesh->n_faces());
//	mesh->write_file("tmp/mesh.off");

	delete [] V;
	delete [] F;

	// ------------------------------------------------------

	cudaDeviceReset();
/*
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
*/
	TIC(time_ptp)
	// BEGIN PTP

	CHE * h_mesh = new CHE(mesh);
	CHE * dd_mesh, * d_mesh;
	cuda_create_CHE(h_mesh, dd_mesh, d_mesh);

	distance_t * h_dist = new distance_t[h_mesh->n_vertices];

	distance_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(distance_t) * h_mesh->n_vertices);
	cudaMalloc(&d_dist[1], sizeof(distance_t) * h_mesh->n_vertices);

	distance_t * error = iter_error_run_ptp_coalescence_gpu(d_mesh, h_mesh->n_vertices, h_dist, d_dist, sources, limits, inv, exact_dist);
	
	delete [] h_dist;
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cuda_free_CHE(dd_mesh, d_mesh);

	// END PTP
	TOC(time_ptp)
/*
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time_ptp, start, stop);
	time_ptp /= 1000;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
*/

	delete mesh;
	delete [] inv;

	return error;
}

/// Return an array of time in seconds.
double * times_farthest_point_sampling_ptp_coalescence_gpu(che * mesh, vector<index_t> & samples, size_t n, distance_t radio)
{
	cudaDeviceReset();

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// BEGIN FPS PTP
	
	vertex * V = new vertex[mesh->n_vertices()];
	index_t * F = new index_t[mesh->n_faces() * che::P];
	index_t * inv = new index_t[mesh->n_vertices()];


	distance_t * h_dist = new distance_t[mesh->n_vertices()];

	distance_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(distance_t) * mesh->n_vertices());
	cudaMalloc(&d_dist[1], sizeof(distance_t) * mesh->n_vertices());

	vector<index_t> limits;
	index_t * toplesets = new index_t[mesh->n_vertices()];
	index_t * sorted_index = new index_t[mesh->n_vertices()];

	cublasHandle_t handle;
	cublasCreate(&handle);

	if(n >= mesh->n_vertices()) n = mesh->n_vertices() >> 1;

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
		
		// sort data by levels, must be improve the coalescence
	
		#pragma omp parallel for
		for(index_t i = 0; i < mesh->n_vertices(); i++)
		{
			V[i] = mesh->gt(sorted_index[i]);
			inv[sorted_index[i]] = i;
		}

		#pragma omp parallel for
		for(index_t he = 0; he < mesh->n_half_edges(); he++)
			F[he] = inv[mesh->vt(he)];

		che * tmp_mesh = new che_off(V, mesh->n_vertices(), F, mesh->n_faces());

		CHE * h_mesh = new CHE(tmp_mesh);
		CHE * dd_mesh, * d_mesh;
		cuda_create_CHE(h_mesh, dd_mesh, d_mesh);

		// exec algorithm
		d = run_ptp_coalescence_gpu(d_mesh, h_mesh->n_vertices, h_dist, d_dist, samples, limits, inv);

		// free memory
		cuda_free_CHE(dd_mesh, d_mesh);
		delete tmp_mesh;

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
		
		samples.push_back(sorted_index[f - 1]);
	}

	cublasDestroy(handle);
	
	delete [] V;
	delete [] F;
	delete [] inv;
	delete [] h_dist;
	delete [] toplesets;
	delete [] sorted_index;

	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);

	// END FPS PTP

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return times;
}

distance_t * iter_error_run_ptp_coalescence_gpu(CHE * d_mesh, const index_t & n_vertices, distance_t * h_dist, distance_t ** d_dist, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * inv, const distance_t * exact_dist)
{
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; v++)
		h_dist[v] = INFINITY;

	for(index_t i = 0; i < sources.size(); i++)
		h_dist[inv[sources[i]]] = 0;

	cudaMemcpy(d_dist[0], h_dist, sizeof(distance_t) * n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dist[1], h_dist, sizeof(distance_t) * n_vertices, cudaMemcpyHostToDevice);

	index_t d = 0, e = 0;
	index_t start, end;
	index_t iter = iterations(limits);

	distance_t * dist_error = new distance_t[iter - limits.size()];	

	for(index_t i = 2; i < iter; i++)
	{
		start = start_v(i, limits);
		end = end_v(i, limits);

		relax_ptp_coalescence <<< NB(end - start), NT >>> (d_mesh, d_dist[!d], d_dist[d], end, start);
		cudaMemcpy(h_dist, d_dist[!d], sizeof(distance_t) * n_vertices, cudaMemcpyDeviceToHost);
		
		// calculating iteration error
		if(i >= limits.size())
		{
			distance_t & error = dist_error[e++] = 0;

			#pragma omp parallel for reduction(+: error)
			for(index_t v = 0; v < n_vertices; v++)
				if(exact_dist[v] > 0)
					error += abs(h_dist[inv[v]] - exact_dist[v]) / exact_dist[v];

			error /= n_vertices - sources.size();
		}

		d = !d;
	}

	return dist_error;
}

