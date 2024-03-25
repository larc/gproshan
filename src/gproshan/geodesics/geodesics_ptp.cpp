#include <gproshan/geodesics/geodesics_ptp.h>

#include <cmath>
#include <cstring>
#include <cassert>


// geometry processing and shape analysis framework
namespace gproshan {


ptp_out_t::ptp_out_t(float *const d, index_t *const c): dist(d), clusters(c) {}


void parallel_toplesets_propagation_cpu(const ptp_out_t & ptp_out,
										const che * mesh,
										const std::vector<index_t> & sources,
										const toplesets_t & toplesets,
										const bool coalescence,
										const bool set_inf,
										const f_ptp<float> & fun
										)
{
	const size_t n_vertices = mesh->n_vertices;

	che * h_mesh = nullptr;
	index_t * inv = nullptr;
	if(coalescence)
	{
		h_mesh = new che(*mesh, toplesets.index, {false, false, false});
		inv = new index_t[n_vertices];

		#pragma omp parallel for
		for(index_t i = 0; i < toplesets.limits.back(); ++i)
			inv[toplesets.index[i]] = i;
	}


	float * dist[2] = {	coalescence ? new float[n_vertices] : ptp_out.dist,
							new float[n_vertices]
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

	const index_t i = run_ptp(	h_mesh ? h_mesh : mesh, sources, toplesets.limits, dist, clusters,
								coalescence ? inv : toplesets.index,
								coalescence ? nullptr : (index_t *) toplesets.index,
								fun);

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		dist[!i][v] = dist[i][v];

	if(coalescence)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
			ptp_out.dist[v] = dist[0][inv[v]];

		delete [] dist[0];
	}

	delete [] dist[1];
	delete [] clusters[1];

	delete [] inv;
	delete h_mesh;
}

void normalize_ptp(float * dist, const size_t n)
{
	float max_d = 0;

	#pragma omp parallel for reduction(max: max_d)
	for(index_t v = 0; v < n; ++v)
		if(dist[v] < INFINITY)
			max_d = std::max(dist[v], max_d);

	#pragma omp parallel for
	for(index_t v = 0; v < n; ++v)
		dist[v] /= max_d;
}


} // namespace gproshan

