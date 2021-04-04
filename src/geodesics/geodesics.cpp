#include "geodesics/geodesics.h"
#include "geodesics/geodesics_ptp.h"

#include "geodesics/heat_method.h"

#include <queue>
#include <cassert>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


geodesics::geodesics(che * mesh, const vector<index_t> & sources, const params & p): n_vertices(mesh->n_vertices)
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

	assert(sources.size() > 0);
	execute(mesh, sources, p);
}

geodesics::~geodesics()
{
	if(free_dist) delete [] dist;

	delete [] sorted_index;
	delete [] clusters;
}

const real_t & geodesics::operator[](const index_t & i) const
{

	assert(i < n_vertices);
	return dist[i];
}

const index_t & geodesics::operator()(const index_t & i) const
{
	assert(i < n_vertices);
	return sorted_index[i];
}

const real_t & geodesics::radio() const
{
	assert(n_sorted != 0);
	return dist[farthest()];
}

const index_t & geodesics::farthest() const
{
	assert(n_sorted != 0);
	return sorted_index[n_sorted - 1];
}

const size_t & geodesics::n_sorted_index() const
{
	return n_sorted;
}

void geodesics::copy_sorted_index(index_t * indexes, const size_t & n) const
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
		dist[sorted_index[i]] /= max;
}

void geodesics::execute(che * mesh, const vector<index_t> & sources, const params & p)
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

void geodesics::run_fastmarching(che * mesh, const vector<index_t> & sources, const size_t & n_iter, const real_t & radio, const fm_function_t & fun)
{
	index_t BLACK = 0, GREEN = 1, RED = 2;
	index_t * color = new index_t[n_vertices];

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		color[v] = GREEN;

	size_t green_count = n_iter ? n_iter : n_vertices;

	priority_queue<pair<real_t, size_t>,
			vector<pair<real_t, size_t> >,
			greater<pair<real_t, size_t> > > Q;

	real_t dv, dp;
	index_t dir; // dir propagation
	vertex vx;

	size_t black_i, v;

	index_t c = 0;
	n_sorted = 0;
	for(index_t s: sources)
	{
		dist[s] = 0;
		if(clusters) clusters[s] = ++c;
		color[s] = RED;
		Q.push(make_pair(dist[s], s));
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

		link_t black_link;
		mesh->link(black_link, black_i);
		for(const index_t & he: black_link)
		{
			v = mesh->vt(he);

			if(color[v] == GREEN)
				color[v] = RED;

			if(color[v] == RED)
			{
				dv = dist[v];
				for_star(v_he, mesh, v)
				{
					dp = update(dir, mesh, v_he, vx);
					//dp = update_step(mesh, dist, v_he);
					if(dp < dv)
					{
						dv = dp;

						if(clusters)
							clusters[v] = dist[mesh->vt(prev(v_he))] < dist[mesh->vt(next(he))] ? clusters[mesh->vt(prev(he))] : clusters[mesh->vt(next(he))];
					}
				}

				if(dv < dist[v])
					Q.push(make_pair(dist[v] = dv, v));
			}
		}
	}

	delete [] color;
}

void geodesics::run_parallel_toplesets_propagation_cpu(che * mesh, const vector<index_t> & sources)
{
	index_t * toplesets = new index_t[n_vertices];
	vector<index_t> limits;
	mesh->compute_toplesets(toplesets, sorted_index, limits, sources);

	double time_ptp;

	TIC(time_ptp)
		parallel_toplesets_propagation_coalescence_cpu({dist, clusters}, mesh, sources, {limits, sorted_index});
	TOC(time_ptp)

	gproshan_log_var(time_ptp);

	delete [] toplesets;
}

void geodesics::run_heat_method(che * mesh, const vector<index_t> & sources)
{
	double time_total, solve_time;
	TIC(time_total)
	solve_time = heat_method(dist, mesh, sources, HEAT_CHOLMOD);
	TOC(time_total)

	gproshan_log_var(time_total - solve_time);
	gproshan_log_var(solve_time);
}


#ifdef GPROSHAN_CUDA

void geodesics::run_parallel_toplesets_propagation_gpu(che * mesh, const vector<index_t> & sources)
{
	index_t * toplesets = new index_t[n_vertices];
	vector<index_t> limits;
	mesh->compute_toplesets(toplesets, sorted_index, limits, sources);

	double time_ptp;
	if(sources.size() > 1)
		time_ptp = parallel_toplesets_propagation_gpu({dist, clusters}, mesh, sources, {limits, sorted_index});
	else
		time_ptp = parallel_toplesets_propagation_coalescence_gpu({dist, clusters}, mesh, sources, {limits, sorted_index});

	gproshan_log_var(time_ptp);

	delete [] toplesets;
}

void geodesics::run_heat_method_gpu(che * mesh, const vector<index_t> & sources)
{
	double time_total, solve_time;
	TIC(time_total)
	solve_time = heat_method(dist, mesh, sources, HEAT_CUDA);
	TOC(time_total)

	gproshan_log_var(time_total - solve_time);
	gproshan_log_var(solve_time);
}

#endif // GPROSHAN_CUDA


//d = {NIL, 0, 1} cross edge, next, prev
real_t geodesics::update(index_t & d, che * mesh, const index_t & he, vertex & vx)
{
	d = NIL;

	a_mat X(3,2);
	index_t x[3];

	x[0] = mesh->vt(next(he));
	x[1] = mesh->vt(prev(he));
	x[2] = mesh->vt(he);				//update x[2]

	vx = mesh->gt(x[2]);

	vertex v[2];
	v[0] = mesh->gt(x[0]) - vx;
	v[1] = mesh->gt(x[1]) - vx;

	X(0, 0) = v[0][0];
	X(1, 0) = v[0][1];
	X(2, 0) = v[0][2];

	X(0, 1) = v[1][0];
	X(1, 1) = v[1][1];
	X(2, 1) = v[1][2];

	return planar_update(d, X, x, vx);
}

real_t geodesics::planar_update(index_t & d, a_mat & X, index_t * x, vertex & vx)
{
	a_mat ones(2,1);
	ones.ones(2,1);

	a_mat Q;
	if(!inv_sympd(Q, X.t() * X))
		return INFINITY;

	a_mat t(2,1);

	t(0) = dist[x[0]];
	t(1) = dist[x[1]];

	real_t p;
	a_mat delta = ones.t() * Q * t;
	real_t dis = as_scalar(delta * delta - (ones.t() * Q * ones) * (as_scalar(t.t() * Q * t) - 1));

	if(dis >= 0)
	{
		p = delta(0) + sqrt(dis);
		p /= as_scalar(ones.t() * Q * ones);
	}
	else p = INFINITY;

	a_mat n = X * Q * (t - p * ones);
	a_mat cond = Q * X.t() * n;

	a_vec v(3);

	if(t(0) == INFINITY || t(1) == INFINITY || dis < 0 || (cond(0) >= 0 || cond(1) >= 0))
	{
		real_t dp[2];
		dp[0] = dist[x[0]] + norm(X.col(0));
		dp[1] = dist[x[1]] + norm(X.col(1));

		d = dp[1] < dp[0];
		v = X.col(d);
		p = dp[d];
	}
	else
	{
		a_mat A(3,2);
		A.col(0) = -n;
		A.col(1) = X.col(1) - X.col(0);
		a_vec b = -X.col(0);
		a_mat l =	solve(A, b);
		v = l(1) * A.col(1) + X.col(0);
	}

	vx += *((vertex *) v.memptr());

	return p;
}


} // namespace gproshan

