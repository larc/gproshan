#include "geodesics_ptp.h"

#include <cmath>

index_t iterations(const vector<index_t> & limits)
{
	return limits.size() + (limits.size() >> 1);
}

index_t start_v(const index_t & i, const vector<index_t> & limits)
{
	return limits[i >> 1];
}

index_t end_v(const index_t & i, const vector<index_t> & limits)
{
	return i < limits.size() ? limits[i] : limits.back();
}

distance_t * parallel_toplesets_propagation_cpu(che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, index_t * clusters)
{
	distance_t * dist[2] = {new distance_t[mesh->n_vertices()], new distance_t[mesh->n_vertices()]};

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices(); v++)
		dist[0][v] = dist[1][v] = INFINITY;

	for(index_t i = 0; i < sources.size(); i++)
	{
		dist[0][sources[i]] = dist[1][sources[i]] = 0;
		if(clusters) clusters[sources[i]] = i + 1;
	}

	index_t d = 1;
	index_t start, end;
	index_t iter = iterations(limits);
	for(index_t i = 2; i < iter; i++)
	{
		start = start_v(i, limits);
		end = end_v(i, limits);

		#pragma omp parallel for
		for(index_t vi = start; vi < end; vi++)
		{
			const index_t & v = sorted_index[vi];
			dist[!d][v] = dist[d][v];

			distance_t p;
			for_star(he, mesh, v)
			{
				p = update_step(mesh, dist[d], he);
				if(p < dist[!d][v])
				{
					dist[!d][v] = p;
					if(clusters)
						clusters[v] = clusters[mesh->vt(prev(he))] != NIL ? clusters[mesh->vt(prev(he))] : clusters[mesh->vt(next(he))];
				}
			}
		}

		d = !d;
	}

	delete [] dist[d];
	return dist[!d];
}

distance_t update_step(che * mesh, const distance_t * dist, const index_t & he)
{
	index_t x[3];
	x[0] = mesh->vt(next(he));
	x[1] = mesh->vt(prev(he));
	x[2] = mesh->vt(he);

	vertex X[2];
	X[0] = mesh->gt(x[0]) - mesh->gt(x[2]);
	X[1] = mesh->gt(x[1]) - mesh->gt(x[2]);

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

	vertex n(tp[0] * (X[0][0]*Q[0][0] + X[1][0]*Q[1][0]) + tp[1] * (X[0][0]*Q[0][1] + X[1][0]*Q[1][1]),
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

void normalize_ptp(distance_t * dist, const size_t & n)
{
	distance_t max_d = 0;

	#pragma omp parallel for reduction(max: max_d)
	for(index_t v = 0; v < n; v++)
		if(dist[v] < INFINITY)
			max_d = max(dist[v], max_d);

	#pragma omp parallel for
	for(index_t v = 0; v < n; v++)
		dist[v] /= max_d;
}

