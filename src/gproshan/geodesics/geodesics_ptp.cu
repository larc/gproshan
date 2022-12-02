#include <gproshan/geodesics/geodesics_ptp.cuh>
#include <gproshan/geodesics/geodesics_ptp.h>

#include <cstdio>
#include <fstream>
#include <cassert>
#include <cublas_v2.h>

#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>


// geometry processing and shape analysis framework
namespace gproshan {


double parallel_toplesets_propagation_gpu(const ptp_out_t & ptp_out, const che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets, const bool & coalescence, const bool & set_inf)
{
	CHE h_mesh(mesh);
	const size_t & n_vertices = h_mesh.n_vertices;

	index_t * inv = nullptr;
	if(coalescence)
	{
		inv = new index_t[n_vertices];
		h_mesh.GT = new vertex[n_vertices];
		h_mesh.EVT = new index_t[n_vertices];
		h_mesh.VT = new index_t[h_mesh.n_half_edges];

		#pragma omp parallel for
		for(index_t i = 0; i < toplesets.limits.back(); ++i)
		{
			h_mesh.GT[i] = mesh->point(toplesets.index[i]);
			inv[toplesets.index[i]] = i;
		}

		#pragma omp parallel for
		for(index_t he = 0; he < mesh->n_half_edges; ++he)
		{
			const index_t & v = mesh->halfedge(he);
			if(v != NIL)
			{
				h_mesh.VT[he] = inv[v];
				if(mesh->evt(v) == he)
					h_mesh.EVT[inv[v]] = he;
			}
		}
	}

	cudaDeviceReset();

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);


	CHE * dd_mesh, * d_mesh;
	cuda_create_CHE(&h_mesh, dd_mesh, d_mesh);

	real_t * h_dist = coalescence ? new real_t[h_mesh.n_vertices] : ptp_out.dist;

	if(set_inf)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < h_mesh.n_vertices; ++v)
			h_dist[v] = INFINITY;
	}

	real_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(real_t) * h_mesh.n_vertices);
	cudaMalloc(&d_dist[1], sizeof(real_t) * h_mesh.n_vertices);

	index_t * d_sorted = nullptr;
	if(!coalescence)
	{
		cudaMalloc(&d_sorted, sizeof(index_t) * h_mesh.n_vertices);
	}

	real_t * d_error;
	cudaMalloc(&d_error, sizeof(real_t) * h_mesh.n_vertices);


	index_t * h_clusters = coalescence && ptp_out.clusters ? new index_t[h_mesh.n_vertices]
															: ptp_out.clusters;

	index_t * d_clusters[2] = {};

	if(h_clusters)
	{
		cudaMalloc(&d_clusters[0], sizeof(index_t) * h_mesh.n_vertices);
		cudaMalloc(&d_clusters[1], sizeof(index_t) * h_mesh.n_vertices);
	}

	const index_t & d = run_ptp_gpu(d_mesh, sources, h_mesh.n_vertices,
									h_dist, d_dist,
									{toplesets.limits, coalescence ? inv : toplesets.index},
									d_error,
									h_clusters, d_clusters,
									d_sorted);

	cudaMemcpy(h_dist, d_dist[d], sizeof(real_t) * h_mesh.n_vertices, cudaMemcpyDeviceToHost);
	if(coalescence)
	{
		#pragma omp parallel for
		for(index_t i = 0; i < toplesets.limits.back(); ++i)
			ptp_out.dist[toplesets.index[i]] = h_dist[i];

		delete [] h_dist;
	}

	if(h_clusters)
	{
		cudaMemcpy(h_clusters, d_clusters[d], sizeof(index_t) * h_mesh.n_vertices, cudaMemcpyDeviceToHost);

		if(coalescence)
		{
			#pragma omp parallel for
			for(index_t i = 0; i < h_mesh.n_vertices; ++i)
				ptp_out.clusters[toplesets.index[i]] = h_clusters[i];

			delete [] h_clusters;
		}

		cudaFree(d_clusters[0]);
		cudaFree(d_clusters[1]);
	}

	cudaFree(d_error);
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cuda_free_CHE(dd_mesh, d_mesh);

	if(coalescence)
	{
		delete [] h_mesh.GT;
		delete [] h_mesh.VT;
		delete [] h_mesh.EVT;
	}
	else
	{
		cudaFree(d_sorted);
	}

	delete [] inv;

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	float time;
	cudaEventElapsedTime(&time, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return time / 1000;
}

index_t run_ptp_gpu(const CHE * d_mesh, const std::vector<index_t> & sources, const index_t & n_vertices,
					real_t * h_dist, real_t ** d_dist, const toplesets_t & inv, real_t * d_error,
					index_t * h_clusters, index_t ** d_clusters, index_t * d_sorted)
{
	for(index_t i = 0; i < sources.size(); ++i)
	{
		const index_t & s = sources[i];
		const index_t & v = d_sorted ? s: inv.index[s];

		h_dist[v] = 0;

		if(h_clusters)
			h_clusters[v] = i;
	}

	cudaMemcpy(d_dist[0], h_dist, sizeof(real_t) * n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dist[1], h_dist, sizeof(real_t) * n_vertices, cudaMemcpyHostToDevice);

	if(d_sorted)
	{
		cudaMemcpy(d_sorted, inv.index, sizeof(index_t) * n_vertices, cudaMemcpyHostToDevice);
	}


	if(h_clusters)
	{
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

		h_clusters ? relax_ptp<<< NB(end - start), NT >>>(d_mesh, d_dist[!d], d_dist[d], d_clusters[!d], d_clusters[d], start, end, d_sorted)
					: relax_ptp<<< NB(end - start), NT >>>(d_mesh, d_dist[!d], d_dist[d], nullptr, nullptr, start, end, d_sorted);

		cudaDeviceSynchronize();

		relative_error<<< NB(n_cond), NT >>>(d_error, d_dist[!d], d_dist[d], start, start + n_cond);
		cudaDeviceSynchronize();

		if(n_cond == thrust::count_if(thrust::device, d_error + start, d_error + start + n_cond, is_ok()))
			++i;

		if(j < inv.limits.size() - 1) ++j;

		d = !d;
	}

	return d;
}

real_t farthest_point_sampling_ptp_gpu(che * mesh, std::vector<index_t> & samples, double & time_fps, size_t n, real_t radio)
{
	cudaDeviceReset();

	float time;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	// BEGIN FPS PTP

	CHE * h_mesh = new CHE(mesh);
	CHE * dd_mesh, * d_mesh;
	cuda_create_CHE(h_mesh, dd_mesh, d_mesh);

	real_t * h_dist = new real_t[h_mesh->n_vertices];
	#pragma omp parallel for
	for(index_t v = 0; v < h_mesh->n_vertices; ++v)
		h_dist[v] = INFINITY;


	real_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(real_t) * h_mesh->n_vertices);
	cudaMalloc(&d_dist[1], sizeof(real_t) * h_mesh->n_vertices);

	real_t * d_error;
	cudaMalloc(&d_error, sizeof(real_t) * h_mesh->n_vertices);

	index_t * d_sorted;
	cudaMalloc(&d_sorted, sizeof(index_t) * h_mesh->n_vertices);

	std::vector<index_t> limits;
	index_t * toplesets = new index_t[h_mesh->n_vertices];
	index_t * sorted_index = new index_t[h_mesh->n_vertices];

	cublasHandle_t handle;
	cublasCreate(&handle);

	if(n >= h_mesh->n_vertices) n = h_mesh->n_vertices >> 1;

	n -= samples.size();
	samples.reserve(n);

	index_t d;
	int f;
	real_t max_dist = INFINITY;
	while(n-- && max_dist > radio)
	{
		limits.clear();
		mesh->compute_toplesets(toplesets, sorted_index, limits, samples);

		d = run_ptp_gpu(d_mesh, samples, h_mesh->n_vertices, h_dist, d_dist, {limits, sorted_index}, d_error, nullptr, nullptr, d_sorted);

		// 1 indexing
		#ifdef GPROSHAN_FLOAT
			cublasIsamax(handle, mesh->n_vertices, d_dist[d], 1, &f);
		#else
			cublasIdamax(handle, mesh->n_vertices, d_dist[d], 1, &f);
		#endif

		if(radio > 0 || !n)
			cudaMemcpy(&max_dist, d_dist[d] + f - 1, sizeof(real_t), cudaMemcpyDeviceToHost);

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
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	time_fps = time / 1000;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return max_dist;
}

__global__
void relax_ptp(const CHE * mesh, real_t * new_dist, real_t * old_dist, index_t * new_clusters, index_t * old_clusters, const index_t start, const index_t end, const index_t * sorted)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;

	if(v < end)
		relax_ptp(mesh, new_dist, old_dist, new_clusters, old_clusters, sorted ? sorted[v] : v);
}

__global__
void relative_error(real_t * error, const real_t * new_dist, const real_t * old_dist, const index_t start, const index_t end, const index_t * sorted)
{
	index_t i = blockDim.x * blockIdx.x + threadIdx.x + start;

	if(i < end)
	{
		index_t v = sorted ? sorted[i] : i;

		#ifdef GPROSHAN_FLOAT
			error[i] = fabsf(new_dist[v] - old_dist[v]) / old_dist[v];
		#else
			error[i] = fabs(new_dist[v] - old_dist[v]) / old_dist[v];
		#endif
	}
}

__host__ __device__
bool is_ok::operator()(const real_t & val) const
{
	return val < PTP_TOL;
}


} // namespace gproshan

