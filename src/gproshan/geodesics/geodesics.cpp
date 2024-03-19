#include <gproshan/geodesics/geodesics.h>
#include <gproshan/geodesics/geodesics_ptp.h>

#include <gproshan/geodesics/heat_method.h>

#include <queue>
#include <cassert>


// geometry processing and shape analysis framework
namespace gproshan {


geodesics::geodesics(che * mesh, const std::vector<index_t> & sources, const params & p): n_vertices(mesh->n_vertices)
{
	assert(n_vertices > 0);

	free_dist = p.dist_alloc == nullptr;

	dist = free_dist ? new real_t[n_vertices] : p.dist_alloc;
	clusters = p.cluster ? new index_t[n_vertices] : nullptr;
	sorted_index = new index_t[n_vertices];

	n_sorted = 0;

	memset(sorted_index, -1, n_vertices * sizeof(index_t));
	for(index_t v = 0; v < n_vertices; ++v)
		dist[v] = INFINITY;

	assert(size(sources) > 0);
	execute(mesh, sources, p);
}

geodesics::~geodesics()
{
	if(free_dist) delete [] dist;

	delete [] sorted_index;
	delete [] clusters;
}

geodesics::operator const real_t * () const
{
	return dist;
}

real_t geodesics::operator[](const index_t i) const
{
	assert(i < n_vertices);
	return dist[i];
}

index_t geodesics::operator()(const index_t i) const
{
	assert(i < n_vertices);
	return sorted_index[i];
}

real_t geodesics::radio() const
{
	assert(n_sorted != 0);
	return dist[farthest()];
}

index_t geodesics::farthest() const
{
	assert(n_sorted != 0);
	return sorted_index[n_sorted - 1];
}

size_t geodesics::n_sorted_index() const
{
	return n_sorted;
}

void geodesics::copy_sorted_index(index_t * indexes, const size_t n) const
{
	assert(n <= n_sorted);
	memcpy(indexes, sorted_index, n * sizeof(index_t));
}

void geodesics::normalize()
{
	if(!n_sorted)
	{
		normalize_ptp(dist, n_vertices);
		return;
	}

	real_t max = dist[farthest()];

	#pragma omp parallel for
	for(size_t i = 0; i < n_sorted; ++i)
		dist[i] /= max;
}

void geodesics::execute(che * mesh, const std::vector<index_t> & sources, const params & p)
{
	switch(p.alg)
	{
		case FM: run_fastmarching(mesh, sources, p.n_iter, p.radio, p.fun);
			break;
		case PTP_CPU: run_parallel_toplesets_propagation_cpu(mesh, sources);
			break;
		case HEAT_METHOD: run_heat_method(mesh, sources);
			break;

#ifdef GPROSHAN_CUDA
		case PTP_GPU: run_parallel_toplesets_propagation_gpu(mesh, sources);
			break;
		case HEAT_METHOD_GPU: run_heat_method_gpu(mesh, sources);
			break;
#endif // GPROSHAN_CUDA
	}
}

void geodesics::run_fastmarching(che * mesh, const std::vector<index_t> & sources, const size_t n_iter, const real_t radio, const fm_function_t & fun)
{
	CHE cmesh(mesh);

	index_t BLACK = 0, GREEN = 1, RED = 2;
	index_t * color = new index_t[n_vertices];

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		color[v] = GREEN;

	size_t green_count = n_iter ? n_iter : n_vertices;

	std::priority_queue<std::pair<real_t, size_t>,
			std::vector<std::pair<real_t, size_t> >,
			std::greater<std::pair<real_t, size_t> > > Q;

	real_t dv, dp;
	vertex vx;

	size_t black_i;

	index_t c = 0;
	n_sorted = 0;
	for(index_t s: sources)
	{
		dist[s] = 0;
		if(clusters) clusters[s] = ++c;
		color[s] = RED;
		Q.push({dist[s], s});
	}

	while(green_count-- && !Q.empty())
	{
		while(!Q.empty() && color[Q.top().second] == BLACK)
			Q.pop();

		if(Q.empty()) break;

		black_i = Q.top().second;
		color[black_i] = BLACK;
		Q.pop();

		if(dist[black_i] > radio) break;

		sorted_index[n_sorted++] = black_i;

		if(fun && !fun(black_i)) break;

		for(const index_t v: mesh->link(black_i))
		{
			if(color[v] == GREEN)
				color[v] = RED;

			if(color[v] == RED)
			{
				dv = dist[v];
				for(const index_t he: mesh->star(v))
				{
					dp = update_step(&cmesh, dist, {mesh->halfedge(he_next(he)),
													mesh->halfedge(he_prev(he)),
													mesh->halfedge(he)
													});
					if(dp < dv)
					{
						dv = dp;

						if(clusters)
							clusters[v] = clusters[mesh->halfedge(he_prev(he))] ? clusters[mesh->halfedge(he_prev(he))] : clusters[mesh->halfedge(he_next(he))];
					}
				}

				if(dv < dist[v])
					Q.push({dist[v] = dv, v});
			}
		}
	}

	delete [] color;
}

void geodesics::run_parallel_toplesets_propagation_cpu(che * mesh, const std::vector<index_t> & sources)
{
	index_t * toplesets = new index_t[n_vertices];
	std::vector<index_t> limits;
	mesh->compute_toplesets(toplesets, sorted_index, limits, sources);

	double time_ptp;

	TIC(time_ptp)
		parallel_toplesets_propagation_cpu({dist, clusters}, mesh, sources, {limits, sorted_index}, size(sources) == 1);
	TOC(time_ptp)

	gproshan_log_var(time_ptp);

	delete [] toplesets;
}

void geodesics::run_heat_method(che * mesh, const std::vector<index_t> & sources)
{
	double time_total, solve_time;
	TIC(time_total)
	solve_time = heat_method(dist, mesh, sources, HEAT_ARMA);
	TOC(time_total)

	gproshan_log_var(time_total - solve_time);
	gproshan_log_var(solve_time);
}


#ifdef GPROSHAN_CUDA

void geodesics::run_parallel_toplesets_propagation_gpu(che * mesh, const std::vector<index_t> & sources)
{
	index_t * toplesets = new index_t[n_vertices];
	std::vector<index_t> limits;
	mesh->compute_toplesets(toplesets, sorted_index, limits, sources);

	double time_ptp = parallel_toplesets_propagation_gpu({dist, clusters}, mesh, sources, {limits, sorted_index});

	gproshan_log_var(time_ptp);

	delete [] toplesets;
}

void geodesics::run_heat_method_gpu(che * mesh, const std::vector<index_t> & sources)
{
	double time_total, solve_time;
	TIC(time_total)
	solve_time = heat_method(dist, mesh, sources, HEAT_CUDA);
	TOC(time_total)

	gproshan_log_var(time_total - solve_time);
	gproshan_log_var(solve_time);
}

#endif // GPROSHAN_CUDA


} // namespace gproshan

