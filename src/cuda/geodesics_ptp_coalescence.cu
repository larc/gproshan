#include "geodesics_ptp_coalescence.cuh"

#include "geodesics_ptp.cuh"
#include "geodesics_ptp.h"

#include "che_off.h"

#include <cstdio>
#include <fstream>
#include <cassert>
#include <cublas_v2.h>

distance_t * parallel_toplesets_propagation_coalescence_gpu(che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, float & time_ptp, index_t * clusters)
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

	index_t d;
	if(clusters)
	{
		index_t * d_clusters[2] = {NULL, NULL};
		cudaMalloc(&d_clusters[0], sizeof(index_t) * h_mesh->n_vertices);
		cudaMalloc(&d_clusters[1], sizeof(index_t) * h_mesh->n_vertices);

		d = run_ptp_coalescence_gpu(d_mesh, h_mesh->n_vertices, h_dist, d_dist, sources, limits, inv, clusters, d_clusters);
		cudaMemcpy(clusters, d_clusters[d], sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);

		cudaFree(d_clusters[0]);
		cudaFree(d_clusters[1]);
	}
	else
	{
		d = run_ptp_coalescence_gpu(d_mesh, h_mesh->n_vertices, h_dist, d_dist, sources, limits, inv);
	}

	cudaMemcpy(h_dist, d_dist[d], sizeof(distance_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);

	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cuda_free_CHE(dd_mesh, d_mesh);

	// END PTP

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time_ptp, start, stop);
	time_ptp /= 1000;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	// 
	delete mesh;
	delete [] inv;

	distance_t * dist = new distance_t[h_mesh->n_vertices];
	for(index_t i = 0; i < h_mesh->n_vertices; i++)
		dist[sorted_index[i]] = h_dist[i];

	return dist;
}

index_t run_ptp_coalescence_gpu(CHE * d_mesh, const index_t & n_vertices, distance_t * h_dist, distance_t ** d_dist, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * inv, index_t * h_clusters, index_t ** d_clusters)
{
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; v++)
		h_dist[v] = INFINITY;

	for(index_t i = 0; i < sources.size(); i++)
		h_dist[inv[sources[i]]] = 0;

	cudaMemcpy(d_dist[0], h_dist, sizeof(distance_t) * n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dist[1], h_dist, sizeof(distance_t) * n_vertices, cudaMemcpyHostToDevice);

	if(h_clusters) // error, need review
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
			relax_ptp_coalescence <<< NB(end - start), NT >>> (d_mesh, d_dist[!d], d_dist[d], d_clusters[!d], d_clusters[d], end, start);
		}
		else
			relax_ptp_coalescence <<< NB(end - start), NT >>> (d_mesh, d_dist[!d], d_dist[d], end, start);
		cudaDeviceSynchronize();
		d = !d;
	}

	return d;
}

__global__
void relax_ptp_coalescence(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t end, index_t start)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;

	if(v < end)
	{
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
void relax_ptp_coalescence(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t * new_clusters, index_t * old_clusters, index_t end, index_t start)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;

	if(v < end)
	{
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

