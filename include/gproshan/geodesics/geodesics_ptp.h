#ifndef GEODESICS_PTP_H
#define GEODESICS_PTP_H

/**
	Parallel Toplesets Propagation algorithm
	topleset: (Topological Level Set)
*/

#include <gproshan/mesh/che.h>

#include <functional>


#ifdef __CUDACC__
	#include <thrust/count.h>
	#include <thrust/device_vector.h>
	#include <thrust/execution_policy.h>

	#define NT 64
	#define NB(x) (x + NT - 1) / NT
#endif // __CUDACC__

#define PTP_TOL 1e-4


// geometry processing and shape analysis framework
namespace gproshan {


#ifdef __CUDACC__

__global__
void relax_ptp(const che * mesh, float * new_dist, float * old_dist, index_t * new_clusters, index_t * old_clusters, const index_t start, const index_t end, const index_t * sorted = nullptr);

__global__
void relative_error(float * error, const float * new_dist, const float * old_dist, const index_t start, const index_t end, const index_t * sorted = nullptr);

struct is_ok
{
	const float * error = nullptr;

	__host_device__
	bool operator()(const float val) const;

	__host_device__
	bool operator()(const index_t val) const;
};

#endif // __CUDACC__


struct ptp_out_t
{
	float * dist = nullptr;
	index_t * clusters = nullptr;

	ptp_out_t(float *const d, index_t *const c = nullptr);
};

struct toplesets_t
{
	const std::vector<index_t> & limits;
	const index_t * index;
};


template<class T>
using f_ptp = std::function<void(T *, index_t, index_t, index_t, index_t)>;


double parallel_toplesets_propagation_gpu(	const ptp_out_t & ptp_out,
											const che * mesh,
											const std::vector<index_t> & sources,
											const toplesets_t & toplesets,
											const bool coalescence = true,
											const bool set_inf = true,
											const f_ptp<float> & fun = nullptr
											);

void parallel_toplesets_propagation_cpu(	const ptp_out_t & ptp_out,
											const che * mesh,
											const std::vector<index_t> & sources,
											const toplesets_t & toplesets,
											const bool coalescence = true,
											const bool set_inf = true,
											const f_ptp<float> & fun = nullptr
											);


float farthest_point_sampling_ptp_gpu(che * mesh, std::vector<index_t> & samples, double & time_fps, size_t n, float radio = 0);

void normalize_ptp(float * dist, const size_t n);


template<class T>
#ifdef __CUDACC__
__forceinline__
#endif
__host_device__
float update_step(const che * mesh, const T * dist, const uvec3 & x)
{
	const vec<T, 3> X[2] = {mesh->point(x[0]) - mesh->point(x[2]),
							mesh->point(x[1]) - mesh->point(x[2])
							};

	const vec<T, 2> t = {dist[x[0]], dist[x[1]]};

	mat<T, 2> q;
	q[0][0] = dot(X[0], X[0]);
	q[0][1] = dot(X[0], X[1]);
	q[1][0] = dot(X[1], X[0]);
	q[1][1] = dot(X[1], X[1]);

	const T det = q[0][0] * q[1][1] - q[0][1] * q[1][0];

	mat<T, 2> Q;
	Q[0][0] = q[1][1] / det;
	Q[0][1] = -q[0][1] / det;
	Q[1][0] = -q[1][0] / det;
	Q[1][1] = q[0][0] / det;

	const T delta = t[0] * (Q[0][0] + Q[1][0]) + t[1] * (Q[0][1] + Q[1][1]);
	const T dis = delta * delta -
								(Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1]) *
								(t[0] * t[0] * Q[0][0] + t[0] * t[1] * (Q[1][0] + Q[0][1]) + t[1] * t[1] * Q[1][1] - 1);

	T p = (delta + sqrtf(dis)) / (Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1]);

	const vec<T, 2> tp = t - p;
	const vec<T, 3> n = {	tp[0] * (X[0][0]*Q[0][0] + X[1][0]*Q[1][0]) + tp[1] * (X[0][0]*Q[0][1] + X[1][0]*Q[1][1]),
							tp[0] * (X[0][1]*Q[0][0] + X[1][1]*Q[1][0]) + tp[1] * (X[0][1]*Q[0][1] + X[1][1]*Q[1][1]),
			 				tp[0] * (X[0][2]*Q[0][0] + X[1][2]*Q[1][0]) + tp[1] * (X[0][2]*Q[0][1] + X[1][2]*Q[1][1])
							};

	const vec<T, 2> c = Q * vec<T, 2>{dot(X[0], n), dot(X[1], n)};

	if(t[0] == INFINITY || t[1] == INFINITY || dis < 0 || c[0] >= 0 || c[1] >= 0)
	{
		const vec<T, 2> & dp = {dist[x[0]] + norm(X[0]), dist[x[1]] + norm(X[1])};
		p = dp[dp[1] < dp[0]];
	}

	return p;
}


template<class T>
#ifdef __CUDACC__
__forceinline__
#endif
__host_device__
void relax_ptp(const che * mesh, T * new_dist, T * old_dist, index_t * new_clusters, index_t * old_clusters, const index_t v)
{
	float & ndv = new_dist[v] = old_dist[v];
	if(new_clusters) new_clusters[v] = old_clusters[v];

	for(const index_t he: mesh->star(v))
	{
		const uvec3 i = {	mesh->halfedge(he_next(he)),
							mesh->halfedge(he_prev(he)),
							mesh->halfedge(he)
							};

		float d = update_step(mesh, old_dist, i);

		if(d < ndv)
		{
			ndv = d;
			if(new_clusters)
				new_clusters[v] = old_clusters[old_dist[i.y()] < old_dist[i.x()] ? i.y() : i.x()];
		}
	}
}


template<class T>
#ifdef __CUDACC__
index_t run_ptp(const che * mesh, const std::vector<index_t> & sources,
				const std::vector<index_t> & limits, T * error, T ** dist, index_t ** clusters,
				const index_t * idx, index_t * sorted, const f_ptp<T> & fun = nullptr)
#else
index_t run_ptp(const che * mesh, const std::vector<index_t> & sources,
				const std::vector<index_t> & limits, T ** dist, index_t ** clusters,
				const index_t * idx, index_t * sorted, const f_ptp<T> & fun = nullptr)
#endif
{
#ifdef __CUDACC__
	T * h_dist = dist[2];
	index_t * h_clusters = clusters[2];
#endif

	for(index_t i = 0; i < size(sources); ++i)
	{					// !coalescence ?
		const index_t v = sorted ? sources[i] : idx[sources[i]];

	#ifdef __CUDACC__
		h_dist[v] = 0;
		if(h_clusters) h_clusters[v] = i + 1;
	#else
		dist[0][v] = dist[1][v] = 0;
		if(clusters && clusters[0])
			clusters[0][v] = clusters[1][v] = i + 1;
	#endif
	}

#ifdef __CUDACC__
	const size_t n_vertices = limits.back();
	cudaMemcpy(dist[0], h_dist, sizeof(T) * n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(dist[1], h_dist, sizeof(T) * n_vertices, cudaMemcpyHostToDevice);
	if(sorted)
	{
		cudaMemcpy(sorted, idx, sizeof(index_t) * n_vertices, cudaMemcpyHostToDevice);
	}
	if(clusters)
	{
		cudaMemcpy(clusters[0], h_clusters, sizeof(index_t) * n_vertices, cudaMemcpyHostToDevice);
		cudaMemcpy(clusters[1], h_clusters, sizeof(index_t) * n_vertices, cudaMemcpyHostToDevice);
	}
#endif

	const int max_iter = size(limits) << 1;

	int iter = -1;
	index_t count = 0;
	index_t i = 1;
	index_t j = 2;
	while(i < j && ++iter < max_iter)
	{
		if(i < (j >> 1)) i = (j >> 1);	// K/2 limit band size

		const index_t start		= limits[i];
		const index_t end		= limits[j];
		const index_t n_cond	= limits[i + 1] - start;

		T * new_dist = dist[iter & 1];
		T * old_dist = dist[!(iter & 1)];

		index_t * new_cluster = clusters[iter & 1];
		index_t * old_cluster = clusters[!(iter & 1)];

	#ifdef __CUDACC__
		relax_ptp<<< NB(end - start), NT >>>(mesh, new_dist, old_dist, new_cluster, old_cluster, start, end, sorted);
		cudaDeviceSynchronize();

		relative_error<<< NB(n_cond), NT >>>(error, new_dist, old_dist, start, start + n_cond, sorted);
		cudaDeviceSynchronize();

		count = sorted ? thrust::count_if(thrust::device, sorted + start, sorted + start + n_cond, is_ok{error})
						: thrust::count_if(thrust::device, error + start, error + start + n_cond, is_ok{});
	#else
		#pragma omp parallel for
		for(index_t v = start; v < end; ++v)
			relax_ptp(mesh, new_dist, old_dist, new_cluster, old_cluster, sorted ? sorted[v] : v);


		count = 0;
		#pragma omp parallel for
		for(index_t k = start; k < start + n_cond; ++k)
		{
			const index_t v = sorted ? sorted[k] : k;
			if(std::abs(new_dist[v] - old_dist[v]) / old_dist[v] < PTP_TOL)
			{
				#pragma omp atomic
				++count;
			}
		}
	#endif

		if(fun) fun(new_dist, i, j, start, end);

		if(n_cond == count)			++i;
		if(j < size(limits) - 1) 	++j;
	}

	return !(iter & 1);
}


} // namespace gproshan

#endif // GEODESICS_PTP_H

