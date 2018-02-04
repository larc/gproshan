#include "geodesics.h"

#include <queue>
#define DP 5e-2

using namespace arma;

index_t BLACK = 0, GREEN = 1, RED = 2;

geodesics::geodesics(che * mesh, const vector<index_t> & _sources, size_t _n_iter, distance_t _radio,
					 bool _normalized, bool parallel, bool _spherical): sources(_sources)
{
	n_vertices = mesh->n_vertices();
	n_iter = _n_iter ? _n_iter : n_vertices;

	if(n_vertices)
	{
		distances = new distance_t[n_vertices];
		clusters = new index_t[n_vertices];
		paths = new index_t[n_vertices];
		color = new index_t[n_vertices];
		sort_index = new index_t[n_vertices];
	}
	else
	{
		distances = NULL;
		clusters = NULL;
		paths = NULL;
		color = NULL;
		sort_index = NULL;
	}

	n_pesos = 0;

	spherical = _spherical;

	radio = _radio;

	memset(sort_index, 255, n_vertices * sizeof(index_t));
	for(index_t v = 0; v < n_vertices; v++)
		distances[v] = INFINITY;

	run(mesh);

	normalized = _normalized;
	if(normalized) normalize();
}

geodesics::~geodesics()
{
	if(distances)	delete [] distances;
	if(paths)		delete [] paths;
	if(clusters)	delete [] clusters;
	if(color)		delete [] color;
	if(sort_index)	delete [] sort_index;
}

index_t geodesics::operator[](index_t i)
{
	return sort_index[i];
}

index_t geodesics::farthest()
{
	return sort_index[n_pesos - 1];
}

size_t geodesics::get_n_radio()
{
	return n_pesos;
}

void geodesics::get_sort_indexes(index_t * indexes, size_t K)
{
	if(K > n_vertices)
		return;

	memcpy(indexes, sort_index, K * sizeof(index_t));
}

void geodesics::execute(che * mesh)
{
	run(mesh);
	if(normalized) normalize();
}

void geodesics::run(che * mesh)
{
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; v++)
		color[v] = GREEN;

	size_t green_count = n_iter;

	priority_queue<pair<distance_t, size_t>,
			vector<pair<distance_t, size_t> >,
			greater<pair<distance_t, size_t> > > cola;

	distance_t p;
	index_t d; // dir propagation
	vertex vx;

	size_t black_i, v;

	index_t c = 0;
	n_pesos = 0;
	for(index_t s: sources)
	{
		distances[s] = 0;
		clusters[s] = ++c;
		color[s] = RED;
		cola.push(make_pair(distances[s], s));
	}

//	index_t * updates = new index_t[n_vertices];
//	memset(updates, 0, n_vertices * sizeof(index_t));

	while(green_count-- && !cola.empty())
	{
		while(!cola.empty() && color[cola.top().second] == BLACK)
			cola.pop();

		if(cola.empty()) break;

		black_i = cola.top().second;
		color[black_i] = BLACK;
		cola.pop();

		if(distances[black_i] > radio) break;

		sort_index[n_pesos++] = black_i;

		link_t black_link;
		mesh->link(black_link, black_i);
		for(index_t he: black_link)
		{
			v = mesh->vt(he);

			if(color[v] == GREEN)
				color[v] = RED;

			if(color[v] == RED)
			{
//				updates[v]++;
				for_star(v_he, mesh, v)
				{
					//p = update(d, mesh, v_he, vx);
					p = update_step(mesh, distances, v_he);
					if(p < distances[v])
					{
						distances[v] = p;
						clusters[v] = distances[mesh->vt(prev(v_he))] < distances[mesh->vt(next(he))] ? clusters[mesh->vt(prev(he))] : clusters[mesh->vt(next(he))];

						if(d == 0) paths[v] = next(v_he);
						else if(d == 1) paths[v] = prev(v_he);
						else paths[v] = v_he;
					}
				}
				if(distances[v] < INFINITY)
					cola.push(make_pair(distances[v], v));
			}
		}
	}

/*
	
	index_t * rings = new index_t[n_vertices];
	index_t * sorted = new index_t[n_vertices];
	vector<index_t> limites;

	mesh->sort_by_rings(rings, sorted, limites, sources);

	for(index_t i = 0; i < n_vertices; i++)
		cout << i << " " << rings[i] << " " << updates[i] << " " << (rings[i] + 1)/ 2 << endl;
	
	debug(limites.size())
	delete [] rings;
	delete [] sorted;
*/
}

int max_iter = 50;
void geodesics::path_to(vector<vertex> & v_path, che * mesh, const index_t & v)
{
	if(!max_iter--) return;
	debug_me(paths)
	debug(v)
	v_path.push_back(mesh->gt(v));
	
	for(auto s: sources)
		if(s == v) return;

	index_t he = paths[v];
	if(mesh->vt(he) != v)
	{
		path_to(v_path, mesh, mesh->vt(he));
		return;
	}
	
	index_t d;
	vertex vx;
	update(d, mesh, he, vx);
	path_to(v_path, mesh, vx, mesh->ot(next(he)));
}

// vertex in edge he
void geodesics::path_to(vector<vertex> & v_path, che * mesh, const vertex & v, const index_t & he)
{
	debug_me(paths)
	debug(v)
	v_path.push_back(v);

	index_t x[3];
	x[0] = mesh->vt(next(he));
	x[1] = mesh->vt(prev(he));
	x[2] = mesh->vt(he);

	vertex vt[3];
	vt[0] = mesh->gt(x[0]) - v;
	vt[1] = mesh->gt(x[1]) - v;
	vt[2] = mesh->gt(x[2]) - v;

	distance_t p[2];
	index_t d[2];
	vertex vx[2];
	mat X(3, 2);

	X(0, 0) = vt[0][0];
	X(1, 0) = vt[0][1];
	X(2, 0) = vt[0][2];

	X(0, 1) = vt[1][0];
	X(1, 1) = vt[1][1];
	X(2, 1) = vt[1][2];
	p[0] = planar_update(d[0], X, x, vx[0]);
	
	X(0, 0) = vt[1][0];
	X(1, 0) = vt[1][1];
	X(2, 0) = vt[1][2];

	X(0, 1) = vt[2][0];
	X(1, 1) = vt[2][1];
	X(2, 1) = vt[2][2];
	p[1] = planar_update(d[1], X, x + 1, vx[1]);

	index_t pi = p[1] < p[0];

	if(d[pi] != NIL) path_to(v_path, mesh, x[d[pi] + pi]);
	else path_to(v_path, mesh, vx[pi], mesh->ot(pi ? prev(he) : next(he)));
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

void geodesics::normalize()
{
	distance_t max = distances[farthest()];
	
	#pragma omp parallel for
	for(size_t i = 0; i < n_pesos; i++)
		distances[sort_index[i]] /= max;
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

distance_t * fast_geodesics(che * mesh, index_t * source, length_t n_sources, const vector<index_t> & limites, const index_t * sorted)
{
	length_t n_vertices = mesh->n_vertices();

	distance_t * dist[2];
	dist[0] = new distance_t[n_vertices];
	dist[1] = new distance_t[n_vertices];
	
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; v++)
		dist[0][v] = dist[1][v] = INFINITY;

	for(index_t s = 0; s < n_sources; s++)
		dist[0][source[s]] = dist[1][source[s]] = 0;
	
	index_t v, k, d = 1;
	index_t iter = limites.size();
	distance_t p;
	distance_t ratio;

	for(index_t i = 2; i < iter; i++)
	{
		ratio = (limites[i] - limites[i - 1]) / (limites[i - 1] - limites[i - 2]);
	//	if(ratio < 1) ratio = 1.0 / ratio;
		k = log(1 + ratio) / log(2) + 2;
		while(k--)
		{
			#pragma parallel for private(v, p)
			for(index_t r = limites[i - 1]; r < limites[i]; r++)
			{
				v = sorted[r];
				dist[!d][v] = dist[d][v];

				for_star(he, mesh, v)
				{
					p = update_step(mesh, dist[d], he);
					if(p < dist[!d][v]) dist[!d][v] = p;
				}
			}

			d = !d;
		}
	}

	delete [] dist[d];
	return dist[!d];
}

