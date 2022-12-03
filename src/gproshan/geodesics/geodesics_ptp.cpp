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
			clusters[0][v] = clusters[1][v] = i + 1;
	}

	const int & max_iter = toplesets.limits.size() << 1;

	int iter = -1;
	index_t i = 1;
	index_t j = 2;
	while(i < j && ++iter < max_iter)
	{
		if(i < (j >> 1)) i = (j >> 1); // K/2 limit band size

		const index_t & start	= toplesets.limits[i];
		const index_t & end		= toplesets.limits[j];
		const index_t & n_cond	= toplesets.limits[i + 1] - start;

		real_t *& new_dist = dist[iter & 1];
		real_t *& old_dist = dist[!(iter & 1)];

		index_t *& new_cluster = clusters[iter & 1];
		index_t *& old_cluster = clusters[!(iter & 1)];

		#pragma omp parallel for
		for(index_t v = start; v < end; ++v)
			relax_ptp(&h_mesh, new_dist, old_dist, new_cluster, old_cluster, inv ? v : toplesets.index[v]);

		#pragma omp parallel for
		for(index_t v = start; v < start + n_cond; ++v)
			error[v] = abs(new_dist[v] - old_dist[v]) / old_dist[v];

		index_t count = 0;
		#pragma omp parallel for reduction(+: count)
		for(index_t v = start; v < start + n_cond; ++v)
			count += error[v] < PTP_TOL;

		if(n_cond == count) ++i;
		if(j < toplesets.limits.size() - 1) ++j;
	}

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		dist[iter & 1][v] = dist[!(iter & 1)][v];
/*
	for(index_t v = 0; v < n_vertices; ++v)
		gproshan_error_var(dist[d][v]);
	gproshan_error_var(ptp_out.dist);
	gproshan_error_var(dist[0]);
	*/
	if(inv)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
			ptp_out.dist[v] = dist[1][inv[v]];

		delete [] dist[0];
	}

	delete [] error;
	delete [] dist[1];
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

