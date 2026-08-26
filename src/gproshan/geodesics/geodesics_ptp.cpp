#include <gproshan/geodesics/geodesics_ptp.h>

#include <cmath>
#include <cstring>
#include <cassert>


// geometry processing and shape analysis framework
namespace gproshan {


ptp_out_t::ptp_out_t(float *const d, index_t *const c): dist(d), clusters(c) {}


coalescence_ptp::coalescence_ptp(const che * m, const toplesets & tps)
{
	if(m) mesh = new che(*m, tps, {false, false, false});

	inv.assign(tps.n_vertices, NIL);

	#pragma omp parallel for
	for(index_t i = 0; i < std::size(tps); ++i)
		inv[tps.sorted[i]] = i;
}

coalescence_ptp::~coalescence_ptp()
{
	delete mesh;
}

coalescence_ptp::operator const index_t * () const
{
	return inv.data();
}

coalescence_ptp::operator const che * () const
{
	return mesh;
}

const che * coalescence_ptp::operator -> () const
{
	return mesh;
}


double parallel_toplesets_propagation_cpu(	const ptp_out_t & ptp_out,
											const che * mesh,
											const std::vector<index_t> & sources,
											const toplesets & tps,
											const bool coalescence,
											const bool set_inf,
											const f_ptp<float> & fun
											)
{
	double time;
	TIC(time);


	const coalescence_ptp inv(coalescence ? mesh : nullptr, tps);
	const size_t n_vertices = coalescence ? inv->n_vertices : mesh->n_vertices;

	float * dist[2] = {	ptp_out.dist, new float[n_vertices]};
	index_t * clusters[2] = {};
	if(ptp_out.clusters)
	{
		clusters[0] = ptp_out.clusters;
		clusters[1] = new index_t[n_vertices];
	}

	if(set_inf)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
			dist[0][v] = dist[1][v] = INFINITY;
	}

	const index_t i = run_ptp(	coalescence ? inv : mesh, sources, tps.splits, dist, clusters,
								coalescence ? nullptr : tps.sorted,
								coalescence ? nullptr : (const index_t *) inv,
								fun);

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		dist[!i][v] = dist[i][v];

	#pragma omp parallel for
	for(index_t v = 0; v < std::size(tps); ++v)
		ptp_out.dist[tps.sorted[v]] = dist[1][v];

	delete [] dist[1];
	delete [] clusters[1];


	TOC(time);

	return time;
}

double farthest_point_sampling_ptp_cpu(std::vector<index_t> & samples, const che * mesh, size_t n, const float radio)
{
	double time;
	TIC(time);

	const size_t n_vertices = mesh->n_vertices;

	float * dist[2] = {new float[n_vertices], new float[n_vertices]};
	index_t * clusters[2] = {};

	std::vector<index_t> inv(n_vertices, NIL);

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		dist[0][v] = dist[1][v] = INFINITY;

	if(!size(samples)) samples.push_back(0);

	toplesets tps(mesh, samples);

	if(n >= n_vertices) n = n_vertices >> 2;

	n -= size(samples);
	samples.reserve(n);

	int farthest;
	float max_dist = INFINITY;
	while(n-- && radio < max_dist)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < std::size(tps); ++v)
			inv[tps.sorted[v]] = v;

		const index_t i = run_ptp(mesh, samples, tps.splits, dist, clusters, tps.sorted, inv.data());

		farthest = 0;
		#pragma omp parallel for
		for(index_t v = 1; v < n_vertices; ++v)
		{
			if(dist[i][v] < INFINITY)
			{
				#pragma omp critical
				if(dist[i][v] > dist[i][farthest])
					farthest = v;
			}
		}

		if(radio > 0 || !n)
			max_dist = dist[i][farthest];

		samples.push_back(tps.sorted[farthest]);
		tps.reset(mesh, samples);
	}

	delete [] dist[0];
	delete [] dist[1];

	TOC(time);

	return time;
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

