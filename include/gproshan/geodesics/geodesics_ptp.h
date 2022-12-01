#ifndef GEODESICS_PTP_H
#define GEODESICS_PTP_H

/**
	Parallel Toplesets Propagation algorithm
	topleset: (Topological Level Set)
*/

#include <gproshan/mesh/che.h>

#define PTP_TOL 1e-4


// geometry processing and shape analysis framework
namespace gproshan {


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

che * ptp_coalescence(index_t * & inv, const che * mesh, const toplesets_t & toplesets);

double parallel_toplesets_propagation_coalescence_gpu(const ptp_out_t & ptp_out, const che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets, const bool & set_inf = 1);

double parallel_toplesets_propagation_gpu(const ptp_out_t & ptp_out, che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets);

void parallel_toplesets_propagation_coalescence_cpu(const ptp_out_t & ptp_out, che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets);

void parallel_toplesets_propagation_cpu(const ptp_out_t & ptp_out, che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets);

real_t farthest_point_sampling_ptp_gpu(che * mesh, std::vector<index_t> & samples, double & time_fps, size_t n, real_t radio = 0);

void normalize_ptp(real_t * dist, const size_t & n);

template<class T>
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

	T p = (delta + sqrt(dis)) / (Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1]);

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
__host__ __device__
void relax_ptp(const CHE * mesh, T * new_dist, T * old_dist, index_t * new_clusters, index_t * old_clusters, const index_t v)
{
	real_t & ndv = new_dist[v] = old_dist[v];

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


} // namespace gproshan

#endif // GEODESICS_PTP_H

