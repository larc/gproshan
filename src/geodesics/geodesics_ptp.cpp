#include "geodesics/geodesics_ptp.h"

#include <cmath>
#include <cstring>
#include <cassert>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


ptp_out_t::ptp_out_t(real_t *const & d, index_t *const & c): dist(d), clusters(c) {}

che * ptp_coalescence(index_t * & inv, const che * mesh, const toplesets_t & toplesets)
{
	// sort data by levels, must be improve the coalescence
	
	vector<vertex> V(toplesets.limits.back());
	vector<index_t> F;
	F.reserve(mesh->n_half_edges);
	
	inv = !inv ? new index_t[mesh->n_vertices] : inv;
	memset(inv, -1, sizeof(index_t) * mesh->n_vertices);

	#pragma omp parallel for
	for(index_t i = 0; i < toplesets.limits.back(); ++i)
	{
		V[i] = mesh->gt(toplesets.index[i]);
		inv[toplesets.index[i]] = i;
	}

	for(index_t he = 0; he < mesh->n_half_edges; ++he)
		if(inv[mesh->vt(he)] != NIL && inv[mesh->vt(prev(he))] != NIL && inv[mesh->vt(next(he))] != NIL)
			F.push_back(inv[mesh->vt(he)]);
	
	return new che(V.data(), toplesets.limits.back(), F.data(), F.size() / che::mtrig);
}

void parallel_toplesets_propagation_coalescence_cpu(const ptp_out_t & ptp_out, che * mesh, const vector<index_t> & sources, const toplesets_t & toplesets)
{
	const size_t n_vertices = mesh->n_vertices;

	index_t * inv = nullptr;
	mesh = ptp_coalescence(inv, mesh, toplesets);

	// ------------------------------------------------------
	real_t * pdist[2] = {new real_t[mesh->n_vertices], new real_t[mesh->n_vertices]};
	real_t * error = new real_t[mesh->n_vertices];

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
		pdist[0][v] = pdist[1][v] = INFINITY;

	for(index_t i = 0; i < sources.size(); ++i)
	{
		pdist[0][inv[sources[i]]] = pdist[1][inv[sources[i]]] = 0;
		if(ptp_out.clusters) ptp_out.clusters[inv[sources[i]]] = i + 1;
	}

	index_t d = 0;
	index_t start, end, n_cond, count;
	index_t i = 1, j = 2;
	
	// maximum number of iterations
	index_t iter = 0;
	index_t max_iter = toplesets.limits.size() << 1;

	while(i < j && iter++ < max_iter)
	{
		if(i < (j >> 1)) i = (j >> 1); // K/2 limit band size

		start = toplesets.limits[i];
		end = toplesets.limits[j];
		n_cond = toplesets.limits[i + 1] - start;
		
		#pragma omp parallel for
		for(index_t v = start; v < end; ++v)
		{
			pdist[!d][v] = pdist[d][v];

			real_t p;
			for_star(he, mesh, v)
			{
				p = update_step(mesh, pdist[d], he);
				if(p < pdist[!d][v])
				{
					pdist[!d][v] = p;

					if(ptp_out.clusters)
						ptp_out.clusters[v] = ptp_out.clusters[mesh->vt(prev(he))] != NIL ? ptp_out.clusters[mesh->vt(prev(he))] : ptp_out.clusters[mesh->vt(next(he))];
				}
			}
		}

		#pragma omp parallel for
		for(index_t v = start; v < start + n_cond; ++v)
			error[v] = abs(pdist[!d][v] - pdist[d][v]) / pdist[d][v];

		count = 0;
		#pragma omp parallel for reduction(+: count)
		for(index_t v = start; v < start + n_cond; ++v)
			count += error[v] < PTP_TOL;

		if(n_cond == count) i++;
		if(j < toplesets.limits.size() - 1) j++;

		d = !d;
	}
	
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		ptp_out.dist[v] = inv[v] != NIL ? pdist[!d][inv[v]] : INFINITY;
	
	delete [] error;
	delete [] pdist[0];
	delete [] pdist[1];
	delete [] inv;
	delete mesh;
}

void parallel_toplesets_propagation_cpu(const ptp_out_t & ptp_out, che * mesh, const vector<index_t> & sources, const toplesets_t & toplesets)
{
	real_t * pdist[2] = {ptp_out.dist, new real_t[mesh->n_vertices]};
	real_t * error = new real_t[mesh->n_vertices];

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
		pdist[0][v] = pdist[1][v] = INFINITY;

	for(index_t i = 0; i < sources.size(); ++i)
	{
		pdist[0][sources[i]] = pdist[1][sources[i]] = 0;
		if(ptp_out.clusters) ptp_out.clusters[sources[i]] = i + 1;
	}

	index_t d = 0;
	index_t start, end, n_cond, count;
	index_t i = 1, j = 2;
	
	// maximum number of iterations
	index_t iter = 0;
	index_t max_iter = toplesets.limits.size() << 1;

	while(i < j && iter++ < max_iter)
	{
		if(i < (j >> 1)) i = (j >> 1); // K/2 limit band size

		start = toplesets.limits[i];
		end = toplesets.limits[j];
		n_cond = toplesets.limits[i + 1] - start;
		
		#pragma omp parallel for
		for(index_t vi = start; vi < end; ++vi)
		{
			const index_t & v = toplesets.index[vi];
			pdist[!d][v] = pdist[d][v];

			real_t p;
			for_star(he, mesh, v)
			{
				p = update_step(mesh, pdist[d], he);
				if(p < pdist[!d][v])
				{
					pdist[!d][v] = p;

					if(ptp_out.clusters)
						ptp_out.clusters[v] = ptp_out.clusters[mesh->vt(prev(he))] != NIL ? ptp_out.clusters[mesh->vt(prev(he))] : ptp_out.clusters[mesh->vt(next(he))];
				}
			}
		}

		#pragma omp parallel for
		for(index_t vi = start; vi < start + n_cond; ++vi)
		{
			const index_t & v = toplesets.index[vi];
			error[vi] = abs(pdist[!d][v] - pdist[d][v]) / pdist[d][v];
		}

		count = 0;
		#pragma omp parallel for reduction(+: count)
		for(index_t vi = start; vi < start + n_cond; ++vi)
			count += error[vi] < PTP_TOL;

		if(n_cond == count) i++;
		if(j < toplesets.limits.size() - 1) j++;

		d = !d;
	}
	
	delete [] error;
	
	if(ptp_out.dist != pdist[!d])
	{
		memcpy(ptp_out.dist, pdist[!d], mesh->n_vertices * sizeof(real_t));
		delete [] pdist[!d];
	}
	else delete [] pdist[d];
}

real_t update_step(che * mesh, const real_t * dist, const index_t & he)
{
	index_t x[3];
	x[0] = mesh->vt(next(he));
	x[1] = mesh->vt(prev(he));
	x[2] = mesh->vt(he);

	vertex X[2];
	X[0] = mesh->gt(x[0]) - mesh->gt(x[2]);
	X[1] = mesh->gt(x[1]) - mesh->gt(x[2]);

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
	real_t dis = delta * delta -
					(Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1]) * 
					(t[0] * t[0] * Q[0][0] + t[0] * t[1] * (Q[1][0] + Q[0][1]) + t[1] * t[1] * Q[1][1] - 1);

	real_t p = (delta + sqrt(dis)) / (Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1]);

	real_t tp[2];
	tp[0] = t[0] - p;
	tp[1] = t[1] - p;

	vertex n(tp[0] * (X[0][0]*Q[0][0] + X[1][0]*Q[1][0]) + tp[1] * (X[0][0]*Q[0][1] + X[1][0]*Q[1][1]),
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

void normalize_ptp(real_t * dist, const size_t & n)
{
	real_t max_d = 0;

	#pragma omp parallel for reduction(max: max_d)
	for(index_t v = 0; v < n; ++v)
		if(dist[v] < INFINITY)
			max_d = max(dist[v], max_d);

	#pragma omp parallel for
	for(index_t v = 0; v < n; ++v)
		dist[v] /= max_d;
}


} // namespace gproshan

