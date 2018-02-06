#include "geodesics.h"

#include <queue>
#include <cassert>

#define DP 5e-2

using namespace arma;

index_t BLACK = 0, GREEN = 1, RED = 2;

geodesics::geodesics(che * mesh, const vector<index_t> & sources, const size_t & n_iter, const distance_t & radio, const bool & parallel)
{
	n_vertices = mesh->n_vertices();

	if(n_vertices)
	{
		distances = new distance_t[n_vertices];
		clusters = new index_t[n_vertices];
		color = new index_t[n_vertices];
		sorted_index = new index_t[n_vertices];
	}
	else
	{
		distances = NULL;
		clusters = NULL;
		color = NULL;
		sorted_index = NULL;
	}

	n_sorted = 0;

	memset(sorted_index, -1, n_vertices * sizeof(index_t));
	for(index_t v = 0; v < n_vertices; v++)
		distances[v] = INFINITY;
	
	assert(sources.size() > 0);
	execute(mesh, sources, n_iter, radio, parallel);
}

geodesics::~geodesics()
{
	if(distances)		delete [] distances;
	if(clusters)		delete [] clusters;
	if(color)			delete [] color;
	if(sorted_index)	delete [] sorted_index;
}

const index_t & geodesics::operator[](const index_t & i) const
{
	return sorted_index[i];
}

const index_t & geodesics::farthest() const
{
	return sorted_index[n_sorted - 1];
}

const size_t & geodesics::n_sorted_index() const
{
	return n_sorted;
}

void geodesics::copy_sorted_index(index_t * indexes, const size_t & n) const
{
	assert(n < n_sorted);
	memcpy(indexes, sorted_index, n * sizeof(index_t));
}

void geodesics::normalize()
{
	distance_t max = distances[farthest()];
	
	#pragma omp parallel for
	for(size_t i = 0; i < n_sorted; i++)
		distances[sorted_index[i]] /= max;
}

void geodesics::execute(che * mesh, const vector<index_t> & sources, const size_t & n_iter, const distance_t & radio, const bool & parallel)
{
	if(parallel) return;
	else run_fastmarching(mesh, sources, n_iter, radio);
}

void geodesics::run_fastmarching(che * mesh, const vector<index_t> & sources, const size_t & n_iter, const distance_t & radio)
{
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; v++)
		color[v] = GREEN;

	size_t green_count = n_iter ? n_iter : n_vertices;

	priority_queue<pair<distance_t, size_t>,
			vector<pair<distance_t, size_t> >,
			greater<pair<distance_t, size_t> > > cola;

	distance_t p;
	index_t d; // dir propagation
	vertex vx;

	size_t black_i, v;

	index_t c = 0;
	n_sorted = 0;
	for(index_t s: sources)
	{
		distances[s] = 0;
		clusters[s] = ++c;
		color[s] = RED;
		cola.push(make_pair(distances[s], s));
	}

	while(green_count-- && !cola.empty())
	{
		while(!cola.empty() && color[cola.top().second] == BLACK)
			cola.pop();

		if(cola.empty()) break;

		black_i = cola.top().second;
		color[black_i] = BLACK;
		cola.pop();

		if(distances[black_i] > radio) break;

		sorted_index[n_sorted++] = black_i;

		link_t black_link;
		mesh->link(black_link, black_i);
		for(index_t he: black_link)
		{
			v = mesh->vt(he);

			if(color[v] == GREEN)
				color[v] = RED;

			if(color[v] == RED)
			{
				for_star(v_he, mesh, v)
				{
					//p = update(d, mesh, v_he, vx);
					p = update_step(mesh, distances, v_he);
					if(p < distances[v])
					{
						distances[v] = p;
						clusters[v] = distances[mesh->vt(prev(v_he))] < distances[mesh->vt(next(he))] ? clusters[mesh->vt(prev(he))] : clusters[mesh->vt(next(he))];
					}
				}

				if(distances[v] < INFINITY)
					cola.push(make_pair(distances[v], v));
			}
		}
	}
}

//d = {NIL, 0, 1} cross edge, next, prev
distance_t geodesics::update(index_t & d, che * mesh, const index_t & he, vertex & vx)
{
	d = NIL;

	mat X(3,2);
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

distance_t geodesics::planar_update(index_t & d, mat & X, index_t * x, vertex & vx)
{
	mat ones(2,1);
	ones.ones(2,1);

	mat Q;
	if(!inv_sympd(Q, X.t() * X))
		return INFINITY;

	mat t(2,1);

	t(0) = distances[x[0]];
	t(1) = distances[x[1]];

	distance_t p;
	mat delta = ones.t() * Q * t;
	distance_t dis = as_scalar(delta*delta - (ones.t() * Q * ones) * (as_scalar(t.t() * Q * t) - 1));

	if(dis >= 0)
	{
		p = delta(0) + sqrt(dis);
		p /= as_scalar(ones.t() * Q * ones);
	}
	else p = INFINITY;

	mat n = X * Q * (t - p * ones);
	mat cond = Q * X.t() * n;

	vec v(3);

	if(t(0) == INFINITY || t(1) == INFINITY || dis < 0 || (cond(0) >= 0 || cond(1) >= 0))
	{
		distance_t dp[2];
		dp[0] = distances[x[0]] + norm(X.col(0));
		dp[1] = distances[x[1]] + norm(X.col(1));

		d = dp[1] < dp[0];
		v = X.col(d);
		p = dp[d];
	}
	else
	{
		mat A(3,2);
		A.col(0) = -n;
		A.col(1) = X.col(1) - X.col(0);
		vec b = -X.col(0);
		mat l =	solve(A, b);
		v = l(1) * A.col(1) + X.col(0);
	}
	
	vx += *((vertex *) v.memptr());

	return p;
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

