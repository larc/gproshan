#ifndef GEODESICS_PTP_H
#define GEODESICS_PTP_H

/**
	Parallel Toplesets Propagation algorithm
	topleset: (Topological Level Set)
*/

#include <gproshan/mesh/che.h>


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
void relax_ptp(const CHE * mesh, real_t * new_dist, real_t * old_dist, index_t * new_clusters, index_t * old_clusters, const index_t start, const index_t end, const index_t * sorted = nullptr);

__global__
void relative_error(real_t * error, const real_t * new_dist, const real_t * old_dist, const index_t start, const index_t end, const index_t * sorted = nullptr);

struct is_ok
{
	__host__ __device__
	bool operator()(const real_t & val) const;
};

#endif // __CUDACC__


struct ptp_out_t
{
	real_t * dist;
	index_t * clusters;

	ptp_out_t(real_t *const & d, index_t *const & c = nullptr);
};

struct toplesets_t
{
	const std::vector<index_t> & limits;
	const index_t *const & index;
};

double parallel_toplesets_propagation_gpu(const ptp_out_t & ptp_out, const che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets, const bool & coalescence = true, const bool & set_inf = true);

void parallel_toplesets_propagation_cpu(const ptp_out_t & ptp_out, che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets, const bool & coalescence = false, const bool & set_inf = true);

real_t farthest_point_sampling_ptp_gpu(che * mesh, std::vector<index_t> & samples, double & time_fps, size_t n, real_t radio = 0);

void normalize_ptp(real_t * dist, const size_t & n);


template<class T>
#ifdef __CUDACC__
__forceinline__
#endif
__host__ __device__
real_t update_step(const CHE * mesh, const T * dist, const uvec3 & x)
{
	const vec<T, 3> X[2] = {mesh->GT[x[0]] - mesh->GT[x[2]],
							mesh->GT[x[1]] - mesh->GT[x[2]]
							};

	const vec<T, 2> & t = {dist[x[0]], dist[x[1]]};

	mat<T, 2> q;
	q[0][0] = dot(X[0], X[0]);
	q[0][1] = dot(X[0], X[1]);
	q[1][0] = dot(X[1], X[0]);
	q[1][1] = dot(X[1], X[1]);

	const T & det = q[0][0] * q[1][1] - q[0][1] * q[1][0];

	mat<T, 2> Q;
	Q[0][0] = q[1][1] / det;
	Q[0][1] = -q[0][1] / det;
	Q[1][0] = -q[1][0] / det;
	Q[1][1] = q[0][0] / det;

	const T & delta = t[0] * (Q[0][0] + Q[1][0]) + t[1] * (Q[0][1] + Q[1][1]);
	const T & dis = delta * delta -
								(Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1]) *
								(t[0] * t[0] * Q[0][0] + t[0] * t[1] * (Q[1][0] + Q[0][1]) + t[1] * t[1] * Q[1][1] - 1);

#ifdef GPROSHAN_FLOAT
	T p = (delta + sqrtf(dis)) / (Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1]);
#else
	T p = (delta + sqrt(dis)) / (Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1]);
#endif

	const vec<T, 2> & tp = t - p;
	const vec<T, 3> & n = {	tp[0] * (X[0][0]*Q[0][0] + X[1][0]*Q[1][0]) + tp[1] * (X[0][0]*Q[0][1] + X[1][0]*Q[1][1]),
							tp[0] * (X[0][1]*Q[0][0] + X[1][1]*Q[1][0]) + tp[1] * (X[0][1]*Q[0][1] + X[1][1]*Q[1][1]),
			 				tp[0] * (X[0][2]*Q[0][0] + X[1][2]*Q[1][0]) + tp[1] * (X[0][2]*Q[0][1] + X[1][2]*Q[1][1])
							};

	const vec<T, 2> & c = Q * vec<T, 2>{dot(X[0], n), dot(X[1], n)};

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
__host__ __device__
void relax_ptp(const CHE * mesh, T * new_dist, T * old_dist, index_t * new_clusters, index_t * old_clusters, const index_t & v)
{
	real_t & ndv = new_dist[v] = old_dist[v];
	if(new_clusters) new_clusters[v] = old_clusters[v];

	real_t d;
	for_star(he, mesh, v)
	{
		const uvec3 i = {mesh->VT[he_next(he)], mesh->VT[he_prev(he)], mesh->VT[he]};

		d = update_step(mesh, old_dist, i);

		if(d < ndv)
		{
			ndv = d;
			if(new_clusters)
				new_clusters[v] = old_dist[i.y()] < old_dist[i.x()] ? old_clusters[i.y()] : old_clusters[i.x()];
		}
	}
}

template<class T>
#ifdef __CUDACC__
__forceinline__
#else
inline
#endif
index_t run_ptp(const CHE * mesh, const std::vector<index_t> & sources,
				const std::vector<index_t> & limits, T * error, T ** dist, index_t ** clusters,
				const index_t * idx, index_t * sorted)
{
#ifdef __CUDACC__
	T * h_dist = dist[2];
	index_t * h_clusters = clusters[2];
#endif

	for(index_t i = 0; i < sources.size(); ++i)
	{					// !coalescence ?
		const index_t & v = sorted ? sources[i] : idx[sources[i]];

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
	const size_t & n_vertices = limits.back();
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

	const int & max_iter = limits.size() << 1;

	int iter = -1;
	index_t count = 0;
	index_t i = 1;
	index_t j = 2;
	while(i < j && ++iter < max_iter)
	{
		if(i < (j >> 1)) i = (j >> 1);	// K/2 limit band size

		const index_t & start	= limits[i];
		const index_t & end		= limits[j];
		const index_t & n_cond	= limits[i + 1] - start;

		T *& new_dist = dist[iter & 1];
		T *& old_dist = dist[!(iter & 1)];

		index_t *& new_cluster = clusters[iter & 1];
		index_t *& old_cluster = clusters[!(iter & 1)];

	#ifdef __CUDACC__
		relax_ptp<<< NB(end - start), NT >>>(mesh, new_dist, old_dist, new_cluster, old_cluster, start, end, sorted);
		cudaDeviceSynchronize();

		relative_error<<< NB(n_cond), NT >>>(error, new_dist, old_dist, start, start + n_cond);
		cudaDeviceSynchronize();

		count = thrust::count_if(thrust::device, error + start, error + start + n_cond, is_ok());
	#else
		#pragma omp parallel for
		for(index_t v = start; v < end; ++v)
			relax_ptp(mesh, new_dist, old_dist, new_cluster, old_cluster, sorted ? sorted[v] : v);

		#pragma omp parallel for
		for(index_t i = start; i < start + n_cond; ++i)
		{
			const index_t & v = sorted ? sorted[v] : i;
			error[v] = abs(new_dist[v] - old_dist[v]) / old_dist[v];
		}

		count = 0;
		#pragma omp parallel for reduction(+: count)
		for(index_t v = start; v < start + n_cond; ++v)
		{
			const index_t & v = sorted ? sorted[v] : i;
			count += error[v] < PTP_TOL;
		}
	#endif

		if(n_cond == count)			++i;
		if(j < limits.size() - 1) 	++j;
	}

	return !(iter & 1);
}


} // namespace gproshan

#endif // GEODESICS_PTP_H

