#include "geodesics_ptp.h"
#include "geodesics_ptp.cuh"
#include "che.cuh"

#include <fstream>
#include <cublas_v2.h>

// TEST ITER CODE ACCURACY
distance_t * cuda_fastmarching(CHE * h_mesh, CHE * d_mesh, const index_t * source, length_t source_size, const vector<index_t> & limites, const index_t * h_sorted, index_t * h_clusters, distance_t * real_dist)
{
	FILE * f_iter_error;
	if(real_dist) f_iter_error = fopen("iter_error", "w");


	distance_t * h_dist = new distance_t[h_mesh->n_vertices];

	distance_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(distance_t) * h_mesh->n_vertices);
	cudaMalloc(&d_dist[1], sizeof(distance_t) * h_mesh->n_vertices);

	index_t * d_clusters[2] = {NULL, NULL};
	if(h_clusters) cudaMalloc(&d_clusters[0], sizeof(index_t) * h_mesh->n_vertices);
	if(h_clusters) cudaMalloc(&d_clusters[1], sizeof(index_t) * h_mesh->n_vertices);

	index_t * d_sorted;
	cudaMalloc(&d_sorted, sizeof(index_t) * h_mesh->n_vertices);

	for(index_t i = 0; i < h_mesh->n_vertices; i++)
		h_dist[i] = INFINITY;

	for(index_t i = 0; i < source_size; i++)
	{
		h_dist[source[i]] = 0;
		if(h_clusters) h_clusters[source[i]] = i + 1;
	}

	cudaMemcpy(d_dist[0], h_dist, sizeof(distance_t) * h_mesh->n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dist[1], h_dist, sizeof(distance_t) * h_mesh->n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_sorted, h_sorted, sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyHostToDevice);
	if(h_clusters) cudaMemcpy(d_clusters[0], h_clusters, sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyHostToDevice);
	if(h_clusters) cudaMemcpy(d_clusters[1], h_clusters, sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyHostToDevice);

	index_t d = 1;
	index_t start, end;

	index_t iter = real_dist ? limites.size() << 1 : iterations(limites);
	for(index_t i = 2; i < iter; i++)
//	while(i && !stop)
	{
		start = start_v(i, limites);
		end = end_v(i, limites);

		if(start == end) break;

//		fastmarching_relax<<<NB(end - start), NT>>>(d_mesh, d_dist[!d], d_dist[d], end, d_clusters[!d], d_clusters[d], start, d_sorted);
//		relax_ptp <<< NB(end - start), NT >>> (d_mesh, d_dist[!d], d_dist[d], d_clusters[!d], d_clusters[d], d_sorted, end, start);

		if(real_dist && i >= limites.size())
		{
			cudaMemcpy(h_dist, d_dist[!d], sizeof(distance_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);
			distance_t error = 0;

			#pragma omp parallel for reduction(+: error)
			for(index_t v = 1; v < h_mesh->n_vertices; v++)
				error += abs(h_dist[v] - real_dist[v]) / real_dist[v];

			error /= h_mesh->n_vertices - source_size;
			if(error < INFINITY) fprintf(f_iter_error, "%d %.10f\n", i, error);
		}
		else
			cudaDeviceSynchronize();

		d = !d;
	}

	if(real_dist) fclose(f_iter_error);

	cudaMemcpy(h_dist, d_dist[d], sizeof(distance_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);
	if(h_clusters) cudaMemcpy(h_clusters, d_clusters[d], sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cudaFree(d_sorted);
	if(h_clusters)
	{
		cudaFree(d_clusters[0]);
		cudaFree(d_clusters[1]);
	}
	return h_dist;
}

inline index_t farthest(distance_t * d, size_t n)
{
	index_t f = 0;

	#pragma omp parallel for
	for(index_t v = 0; v < n; v++)
		#pragma omp critical
		if(d[v] < INFINITY && d[f] < d[v])
			f = v;

	return f;
}

distance_t farthest_point_sampling_gpu(vector<index_t> & points, float & time, che * mesh, size_t n, distance_t radio)
{
	debug_me(GEODESICS_PTP)

	CHE * h_mesh;
	CHE * dd_mesh;
	CHE * d_mesh;

	h_mesh = new CHE(mesh);

	cuda_create_CHE(h_mesh, dd_mesh, d_mesh);

	ofstream os((PATH_TEST + "fastmarching/" + mesh->name() + ".fps").c_str());

	size_t n_v = mesh->n_vertices();

	index_t * rings = new index_t[n_v];
	index_t * h_sorted = new index_t[n_v];

	distance_t * h_dist = new distance_t[n_v];

	distance_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(distance_t) * n_v);
	cudaMalloc(&d_dist[1], sizeof(distance_t) * n_v);

	index_t * d_sorted;
	cudaMalloc(&d_sorted, sizeof(index_t) * n_v);

	// ---------------------------------------------------------------------------------------------
	cudaEvent_t start;
	cudaEvent_t stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cublasHandle_t handle;
	cublasCreate(&handle);

	if(n >= mesh->n_vertices())
		n = mesh->n_vertices() / 2;

	n -= points.size();
	points.reserve(n);

	time = 0;
	float time_aux;

	index_t d;
	int f;
	distance_t max_dis = INFINITY;
	while(n-- && max_dis > radio)
	{
		cudaEventRecord(start, 0);

		vector<index_t> limites;
		mesh->compute_toplesets(rings, h_sorted, limites, points);

		d = run_ptp_gpu(d_mesh, n_v, h_dist, d_dist, points, limites, h_sorted, d_sorted);

		// 1 indexing T_T
		#ifdef SINGLE_P
		cublasIsamax(handle, mesh->n_vertices(), d_dist[d], 1, &f);
		#else
		cublasIdamax(handle, mesh->n_vertices(), d_dist[d], 1, &f);
		#endif

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time_aux, start, stop);

		time_aux /= 1000;
		time += time_aux;

		os << points.size() << " " << time_aux << endl;

		if(radio > 0 || !n)
			cudaMemcpy(&max_dis, d_dist[d] + f - 1, sizeof(distance_t), cudaMemcpyDeviceToHost);
		points.push_back(f - 1);
	}


	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cublasDestroy(handle);
	// ---------------------------------------------------------------------------------------------

	os.close();

	delete [] rings;
	delete [] h_sorted;
	delete [] h_dist;

	cuda_free_CHE(dd_mesh, d_mesh);
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cudaFree(d_sorted);

	return max_dis;
}

