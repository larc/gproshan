#include "geodesics/geodesics_ptp_coalescence.cuh"

#include "mesh/che_off.h"

#include <cstdio>
#include <fstream>
#include <cassert>
#include <cublas_v2.h>

#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


double parallel_toplesets_propagation_coalescence_gpu(const ptp_out_t & ptp_out, const che * mesh, const vector<index_t> & sources, const toplesets_t & toplesets, const bool & set_inf)
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

index_t run_ptp_coalescence_gpu(CHE * d_mesh, const index_t & n_vertices, real_t * h_dist, real_t ** d_dist, const vector<index_t> & sources, const toplesets_t & inv, real_t * d_error, index_t * h_clusters, index_t ** d_clusters)
{
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		h_dist[v] = INFINITY;

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

		if(h_clusters)
			relax_ptp_coalescence <<< NB(end - start), NT >>> (d_mesh, d_dist[!d], d_dist[d], d_clusters[!d], d_clusters[d], end, start);
		else
			relax_ptp_coalescence <<< NB(end - start), NT >>> (d_mesh, d_dist[!d], d_dist[d], end, start);

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

__global__
void relax_ptp_coalescence(CHE * mesh, real_t * new_dist, real_t * old_dist, index_t end, index_t start)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;

	if(v < end)
	{
		if(v < mesh->n_vertices)
		{
			new_dist[v] = old_dist[v];

			real_t d;
			cu_for_star(he, mesh, v)
			{
				d = cu_update_step(mesh, old_dist, he);
				if(d < new_dist[v]) new_dist[v] = d;
			}
		}
	}
}


__global__
void relax_ptp_coalescence(CHE * mesh, real_t * new_dist, real_t * old_dist, index_t * new_clusters, index_t * old_clusters, index_t end, index_t start)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;

	if(v < end)
	{
		if(v < mesh->n_vertices)
		{
			new_dist[v] = old_dist[v];
			new_clusters[v] = old_clusters[v];

			real_t d;
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

__forceinline__ __device__
real_t cu_update_step(CHE * mesh, const real_t * dist, const index_t & he)
{
	index_t x[3];
	x[0] = mesh->VT[cu_next(he)];
	x[1] = mesh->VT[cu_prev(he)];
	x[2] = mesh->VT[he];

	vertex_cu X[2];
	X[0] = mesh->GT[x[0]] - mesh->GT[x[2]];
	X[1] = mesh->GT[x[1]] - mesh->GT[x[2]];

	real_t t[2];
	t[0] = dist[x[0]];
	t[1] = dist[x[1]];

	real_t q[2][2];
	q[0][0] = (X[0], X[0]);
	q[0][1] = (X[0], X[1]);
	q[1][0] = (X[1], X[0]);
	q[1][1] = (X[1], X[1]);

	real_t det = q[0][0] * q[1][1] - q[0][1] * q[1][0];
	real_t Q[2][2];
	Q[0][0] = q[1][1] / det;
	Q[0][1] = -q[0][1] / det;
	Q[1][0] = -q[1][0] / det;
	Q[1][1] = q[0][0] / det;

	real_t delta = t[0] * (Q[0][0] + Q[1][0]) + t[1] * (Q[0][1] + Q[1][1]);
	real_t dis = delta * delta - (Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1]) * (t[0]*t[0]*Q[0][0] + t[0]*t[1]*(Q[1][0] + Q[0][1]) + t[1]*t[1]*Q[1][1] - 1);

	real_t p;

	if(dis >= 0)
	{
		#ifdef SINGLE_P
			p = delta + sqrtf(dis);
		#else
			p = delta + sqrt(dis);
		#endif
		p /= Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1];
	}

	real_t tp[2];
	tp[0] = t[0] - p;
	tp[1] = t[1] - p;

	vertex_cu n(tp[0] * (X[0][0]*Q[0][0] + X[1][0]*Q[1][0]) + tp[1] * (X[0][0]*Q[0][1] + X[1][0]*Q[1][1]),
			 tp[0] * (X[0][1]*Q[0][0] + X[1][1]*Q[1][0]) + tp[1] * (X[0][1]*Q[0][1] + X[1][1]*Q[1][1]),
			 tp[0] * (X[0][2]*Q[0][0] + X[1][2]*Q[1][0]) + tp[1] * (X[0][2]*Q[0][1] + X[1][2]*Q[1][1]) );

	real_t cond[2];
	cond[0] = (X[0] , n);
	cond[1] = (X[1] , n);

	real_t c[2];
	c[0] = cond[0] * Q[0][0] + cond[1] * Q[0][1];
	c[1] = cond[0] * Q[1][0] + cond[1] * Q[1][1];

	if(t[0] == INFINITY || t[1] == INFINITY || dis < 0 || c[0] >= 0 || c[1] >= 0)
	{
		real_t dp[2];
		dp[0] = dist[x[0]] + *X[0];
		dp[1] = dist[x[1]] + *X[1];

		p = dp[dp[1] < dp[0]];
	}

	return p;
}


} // namespace gproshan

