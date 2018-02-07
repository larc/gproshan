#include "geodesics_ptp.h"
#include "che.cuh"

#include <cstdio>
#include <fstream>
#include <cassert>
#include <cublas_v2.h>

#define NT 32
#define NB(x) (x + NT - 1) / NT

length_t CHE::band = 6;

__device__
distance_t cu_update_step(CHE * mesh, const distance_t * dist, const index_t & he);

__global__
void relax_ptp(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t * new_clusters, index_t * old_clusters, index_t * sorted, index_t end, index_t start = 0);

__global__
void relax_ptp(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t * sorted, index_t end, index_t start = 0);

index_t run_ptp_gpu(CHE * d_mesh, const index_t & n_vertices, distance_t * h_dist, distance_t ** d_dist, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * h_sorted, index_t * d_sorted, index_t * h_clusters = NULL, index_t ** d_clusters = NULL);

distance_t * parallel_toplesets_propagation_gpu(che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, float & time_ptp, index_t * clusters)
{
	debug_me(GEODESICS_PTP)

	cudaDeviceReset();

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
	
	index_t d;
	if(clusters)
	{
		index_t * d_clusters[2] = {NULL, NULL};
		cudaMalloc(&d_clusters[0], sizeof(index_t) * h_mesh->n_vertices);	
		cudaMalloc(&d_clusters[1], sizeof(index_t) * h_mesh->n_vertices);	
	
		d = run_ptp_gpu(d_mesh, h_mesh->n_vertices, h_dist, d_dist, sources, limits, sorted_index, d_sorted, clusters, d_clusters);
		cudaMemcpy(clusters, d_clusters[d], sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);

		cudaFree(d_clusters[0]);
		cudaFree(d_clusters[1]);
	}
	else
	{
		d = run_ptp_gpu(d_mesh, h_mesh->n_vertices, h_dist, d_dist, sources, limits, sorted_index, d_sorted);
	}

	cudaMemcpy(h_dist, d_dist[d], sizeof(distance_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);
	
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cudaFree(d_sorted);
	cuda_free_CHE(dd_mesh, d_mesh);
	
	// END PTP

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time_ptp, start, stop);
	time_ptp /= 1000;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return h_dist;
}

index_t run_ptp_gpu(CHE * d_mesh, const index_t & n_vertices, distance_t * h_dist, distance_t ** d_dist, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * h_sorted, index_t * d_sorted, index_t * h_clusters, index_t ** d_clusters)
{
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; v++)
		h_dist[v] = INFINITY;

	for(index_t i = 0; i < sources.size(); i++)
		h_dist[sources[i]] = 0;
		
	cudaMemcpy(d_dist[0], h_dist, sizeof(distance_t) * n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dist[1], h_dist, sizeof(distance_t) * n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_sorted, h_sorted, sizeof(index_t) * n_vertices, cudaMemcpyHostToDevice);
	
	if(h_clusters)
	{
		assert(d_clusters);

		for(index_t i = 0; i < sources.size(); i++)
			h_clusters[sources[i]] = i + 1;

		cudaMemcpy(d_clusters[0], h_clusters, sizeof(index_t) * n_vertices, cudaMemcpyHostToDevice);
		cudaMemcpy(d_clusters[1], h_clusters, sizeof(index_t) * n_vertices, cudaMemcpyHostToDevice);
	}

	index_t d = 0;
	index_t start, end;
	index_t iter = iterations(limits);
	for(index_t i = 2; i < iter; i++)
	{
		start = start_v(i, limits);
		end = end_v(i, limits);
		
		if(h_clusters)
		{
			relax_ptp <<< NB(end - start), NT >>> (d_mesh, d_dist[!d], d_dist[d], d_clusters[!d], d_clusters[d], d_sorted, end, start);
		}
		else
			relax_ptp <<< NB(end - start), NT >>> (d_mesh, d_dist[!d], d_dist[d], d_sorted, end, start);
		cudaDeviceSynchronize();
		d = !d;
	}
	
	return d;
}

__global__
void relax_ptp(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t * sorted, index_t end, index_t start)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;
	
	if(v < end)
	{
		v = sorted ? sorted[v] : v;
		if(v < mesh->n_vertices)
		{
			new_dist[v] = old_dist[v];

			distance_t d;
			cu_for_star(he, mesh, v)
			{
				d = cu_update_step(mesh, old_dist, he);
				if(d < new_dist[v]) new_dist[v] = d;
			}
		}
	}
}


__global__
void relax_ptp(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t * new_clusters, index_t * old_clusters, index_t * sorted, index_t end, index_t start)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;
	
	if(v < end)
	{
		v = sorted ? sorted[v] : v;
		if(v < mesh->n_vertices)
		{
			new_dist[v] = old_dist[v];
			new_clusters[v] = old_clusters[v];

			distance_t d;
			cu_for_star(he, mesh, v)
			{
				d = cu_update_step(mesh, old_dist, he);
				if(d < new_dist[v])
				{
					new_dist[v] = d;
					new_clusters[v] = old_dist[mesh->VT[cu_prev(he)]] < old_dist[mesh->VT[cu_next(he)]] ? old_clusters[mesh->VT[cu_prev(he)]] : old_clusters[mesh->VT[cu_next(he)]];
				}
			}
		}
	}
}

__device__
distance_t cu_update_step(CHE * mesh, const distance_t * dist, const index_t & he)
{
	index_t x[3];
	x[0] = mesh->VT[cu_next(he)];
	x[1] = mesh->VT[cu_prev(he)];
	x[2] = mesh->VT[he];

	vertex_cu X[2];
	X[0] = mesh->GT[x[0]] - mesh->GT[x[2]];
	X[1] = mesh->GT[x[1]] - mesh->GT[x[2]];

	distance_t t[2];
	t[0] = dist[x[0]];
	t[1] = dist[x[1]];

	distance_t q[2][2];
	q[0][0] = (X[0], X[0]);
	q[0][1] = (X[0], X[1]);
	q[1][0] = (X[1], X[0]);
	q[1][1] = (X[1], X[1]);
	
	distance_t det = q[0][0] * q[1][1] - q[0][1] * q[1][0];
	distance_t Q[2][2];
	Q[0][0] = q[1][1] / det;
	Q[0][1] = -q[0][1] / det;
	Q[1][0] = -q[1][0] / det;
	Q[1][1] = q[0][0] / det;

	distance_t delta = t[0] * (Q[0][0] + Q[1][0]) + t[1] * (Q[0][1] + Q[1][1]);
	distance_t dis = delta * delta - (Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1]) * (t[0]*t[0]*Q[0][0] + t[0]*t[1]*(Q[1][0] + Q[0][1]) + t[1]*t[1]*Q[1][1] - 1);
	
	distance_t p;

	if(dis >= 0)
	{
		p = delta + sqrt(dis);
		p /= Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1];
	}

	distance_t tp[2];
	tp[0] = t[0] - p;
	tp[1] = t[1] - p;

	vertex_cu n(tp[0] * (X[0][0]*Q[0][0] + X[1][0]*Q[1][0]) + tp[1] * (X[0][0]*Q[0][1] + X[1][0]*Q[1][1]),
			 tp[0] * (X[0][1]*Q[0][0] + X[1][1]*Q[1][0]) + tp[1] * (X[0][1]*Q[0][1] + X[1][1]*Q[1][1]),
			 tp[0] * (X[0][2]*Q[0][0] + X[1][2]*Q[1][0]) + tp[1] * (X[0][2]*Q[0][1] + X[1][2]*Q[1][1]) );

	distance_t cond[2];
	cond[0] = (X[0] , n);
	cond[1] = (X[1] , n);

	distance_t c[2];
	c[0] = cond[0] * Q[0][0] + cond[1] * Q[0][1];
	c[1] = cond[0] * Q[1][0] + cond[1] * Q[1][1];

	if(t[0] == INFINITY || t[1] == INFINITY || dis < 0 || c[0] >= 0 || c[1] >= 0)
	{
		distance_t dp[2];
		dp[0] = dist[x[0]] + *X[0];
		dp[1] = dist[x[1]] + *X[1];

		p = dp[dp[1] < dp[0]];
	}

	return p;
}

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
		mesh->sort_by_rings(rings, h_sorted, limites, points);
		
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

