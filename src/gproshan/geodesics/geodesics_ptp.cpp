#include <gproshan/geodesics/geodesics_ptp.h>

#include <cmath>
#include <cstring>
#include <cassert>


// geometry processing and shape analysis framework
namespace gproshan {


ptp_out_t::ptp_out_t(real_t * d, index_t * c): dist(d), clusters(c) {}


void parallel_toplesets_propagation_cpu(const ptp_out_t & ptp_out, che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets, const bool coalescence, const bool set_inf)
{
	CHE h_mesh(mesh);
	const size_t n_vertices = h_mesh.n_vertices;

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
			const index_t v = mesh->halfedge(he);
			if(v != NIL)
			{
				h_mesh.VT[he] = inv[v];
				if(mesh->evt(v) == he)
					h_mesh.EVT[inv[v]] = he;
			}
		}
	}

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

	const index_t i = run_ptp(&h_mesh, sources, toplesets.limits, dist, clusters,
								coalescence ? inv : toplesets.index,
								coalescence ? nullptr : (index_t *) toplesets.index);

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		dist[!i][v] = dist[i][v];

	gproshan_log_var(inv);
	gproshan_log_var(dist[0]);
	gproshan_log_var(dist[1]);
	gproshan_log_var(ptp_out.dist);

	if(coalescence)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
			ptp_out.dist[v] = dist[0][inv[v]];

		delete [] dist[0];
	}

	delete [] dist[1];
	delete [] clusters[1];

	if(coalescence)
	{
		delete [] h_mesh.GT;
		delete [] h_mesh.VT;
		delete [] h_mesh.EVT;
	}

	delete [] inv;
}

void normalize_ptp(real_t * dist, const size_t n)
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

