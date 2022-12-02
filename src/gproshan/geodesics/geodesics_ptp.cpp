#include <gproshan/geodesics/geodesics_ptp.h>

#include <cmath>
#include <cstring>
#include <cassert>


// geometry processing and shape analysis framework
namespace gproshan {


ptp_out_t::ptp_out_t(real_t *const & d, index_t *const & c): dist(d), clusters(c) {}


void parallel_toplesets_propagation_cpu(const ptp_out_t & ptp_out, che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets, const bool & coalescence, const bool & set_inf)
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

	real_t * error = new real_t[n_vertices];
	real_t * dist[2] = {	coalescence ? new real_t[n_vertices] : ptp_out.dist,
							new real_t[n_vertices]
							};
	index_t * clusters[2] = {	coalescence && ptp_out.clusters ? new index_t[n_vertices] : ptp_out.clusters,
								ptp_out.clusters ? new index_t[n_vertices] : nullptr
								};

	if(set_inf)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
			dist[0][v] = dist[1][v] = INFINITY;
	}

	for(index_t i = 0; i < sources.size(); ++i)
	{
		const index_t & s = sources[i];
		const index_t & v = inv ? inv[s] : s;

		dist[0][v] = dist[1][v] = 0;

		if(ptp_out.clusters)
			clusters[0][inv[s]] = clusters[0][inv[s]] = i + 1;
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
			relax_ptp(&h_mesh, dist[!d], dist[d], clusters[!d], clusters[d], inv ? inv[v] : v);

		#pragma omp parallel for
		for(index_t v = start; v < start + n_cond; ++v)
			error[v] = abs(dist[!d][v] - dist[d][v]) / dist[d][v];

		count = 0;
		#pragma omp parallel for reduction(+: count)
		for(index_t v = start; v < start + n_cond; ++v)
			count += error[v] < PTP_TOL;

		if(n_cond == count) ++i;
		if(j < toplesets.limits.size() - 1) ++j;

		d = !d;
	}

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		dist[!d][v] = dist[d][v];

	if(inv)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
			dist[!d][v] = dist[d][inv[v]];
	}

	delete [] error;
	delete [] dist[0];
	delete [] dist[1];
	delete [] clusters[0];
	delete [] clusters[1];
	delete [] inv;
}

void normalize_ptp(real_t * dist, const size_t & n)
{
	real_t max_d = 0;

	#pragma omp parallel for reduction(max: max_d)
	for(index_t v = 0; v < n; ++v)
		if(dist[v] < INFINITY)
			max_d = std::max(dist[v], max_d);

	#pragma omp parallel for
	for(index_t v = 0; v < n; ++v)
		dist[v] /= max_d;
}


} // namespace gproshan

