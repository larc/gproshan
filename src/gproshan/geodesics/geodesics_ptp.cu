#include <gproshan/geodesics/geodesics_ptp.h>

#include <gproshan/mesh/che.cuh>

#include <cstdio>
#include <fstream>
#include <cassert>
#include <cublas_v2.h>


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

	real_t * h_dist = coalescence ? new real_t[n_vertices] : ptp_out.dist;
	index_t * h_clusters = coalescence && ptp_out.clusters ? new index_t[n_vertices]
															: ptp_out.clusters;

	real_t * d_error = nullptr;
	real_t * d_dist[3] = {};
	index_t * d_clusters[3] = {};
	index_t * d_sorted = nullptr;

	cudaMalloc(&d_error, sizeof(real_t) * n_vertices);
	cudaMalloc(&d_dist[0], sizeof(real_t) * n_vertices);
	cudaMalloc(&d_dist[1], sizeof(real_t) * n_vertices);
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

	const index_t & i = run_ptp(d_mesh, sources, toplesets.limits, d_error, d_dist, d_clusters,
								coalescence ? inv : toplesets.index, d_sorted);

	cudaMemcpy(h_dist, d_dist[i], sizeof(real_t) * n_vertices, cudaMemcpyDeviceToHost);

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
	cuda_free_CHE(dd_mesh, d_mesh);

	if(coalescence)
	{
		delete [] h_mesh.GT;
		delete [] h_mesh.VT;
		delete [] h_mesh.EVT;
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

real_t farthest_point_sampling_ptp_gpu(che * mesh, std::vector<index_t> & samples, double & time_fps, size_t n, real_t radio)
{
	CHE h_mesh(mesh);
	const size_t & n_vertices = h_mesh.n_vertices;

	cudaDeviceReset();

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	CHE * dd_mesh, * d_mesh;
	cuda_create_CHE(&h_mesh, dd_mesh, d_mesh);

	real_t * h_dist = new real_t[n_vertices];

	real_t * d_error = nullptr;
	real_t * d_dist[3] = {};
	index_t * d_clusters[3] = {};
	index_t * d_sorted = nullptr;

	cudaMalloc(&d_error, sizeof(real_t) * n_vertices);
	cudaMalloc(&d_dist[0], sizeof(real_t) * n_vertices);
	cudaMalloc(&d_dist[1], sizeof(real_t) * n_vertices);
	cudaMalloc(&d_sorted, sizeof(index_t) * n_vertices);
	d_dist[2] = h_dist;

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		h_dist[v] = INFINITY;

	std::vector<index_t> limits;
	index_t * toplesets = new index_t[n_vertices];
	index_t * sorted_index = new index_t[n_vertices];

	cublasHandle_t handle;
	cublasCreate(&handle);

	if(n >= n_vertices) n = n_vertices >> 2;

	n -= samples.size();
	samples.reserve(n);

	int farthest;
	real_t max_dist = INFINITY;
	while(n-- && radio < max_dist)
	{
		limits.clear();
		mesh->compute_toplesets(toplesets, sorted_index, limits, samples);

		const index_t & i = run_ptp(d_mesh, samples, limits, d_error, d_dist, d_clusters, sorted_index, d_sorted);

		// 1 indexing
		#ifdef GPROSHAN_FLOAT
			cublasIsamax(handle, mesh->n_vertices, d_dist[i], 1, &farthest);
		#else
			cublasIdamax(handle, mesh->n_vertices, d_dist[i], 1, &farthest);
		#endif

		if(radio > 0 || !n)
			cudaMemcpy(&max_dist, d_dist[i] + farthest - 1, sizeof(real_t), cudaMemcpyDeviceToHost);

		samples.push_back(farthest - 1);
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
void relax_ptp(const CHE * mesh, real_t * new_dist, real_t * old_dist, index_t * new_clusters, index_t * old_clusters, const index_t start, const index_t end, const index_t * sorted)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;

	if(v < end)
		relax_ptp(mesh, new_dist, old_dist, new_clusters, old_clusters, sorted ? sorted[v] : v);
}

__global__
void relative_error(real_t * error, const real_t * new_dist, const real_t * old_dist, const index_t start, const index_t end, const index_t * sorted)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;

	if(v < end)
	{
		v = sorted ? sorted[v] : v;

		#ifdef GPROSHAN_FLOAT
			error[v] = fabsf(new_dist[v] - old_dist[v]) / old_dist[v];
		#else
			error[v] = fabs(new_dist[v] - old_dist[v]) / old_dist[v];
		#endif
	}
}

__host_device__
bool is_ok::operator()(const real_t & val) const
{
	return val < PTP_TOL;
}

__host_device__
bool is_ok::operator()(const index_t & i) const
{
	return error[i] < PTP_TOL;
}


} // namespace gproshan

