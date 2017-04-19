#include "fastmarching.h"
#include <queue>

using namespace arma;

size_t black = 0, green = 1, red = 2;

path::path(face fx_, vertex vx_, size_t vi_)
{
	f = fx_;
	v = vx_;
	i = vi_;
}

fastmarching::fastmarching(bool _spherical)
{
	n_vertices = 0;
	distances = 0;
	paths = 0;
	color = 0;
	spherical = _spherical;
	radio = INFINITY;
	sort_index = 0;
	n_pesos = 0;
}

fastmarching::fastmarching(off & model, vector<index_t> & source, distance_t r, bool normalized, bool _spherical)
{
	n_vertices = model.get_nvertices();
	distances = new distance_t[n_vertices];
	cluster = new index_t[n_vertices];
	paths = new path[n_vertices];
	color = new size_t[n_vertices];
	sort_index = new size_t[n_vertices];
	n_pesos = 0;

	spherical = _spherical;

	radio = r;

	calcular_FM(model, source);
	if(normalized) normalize();
}

fastmarching::~fastmarching()
{
	if(distances) delete [] distances;
	if(paths) delete [] paths;
	if(color) delete [] color;
	if(sort_index) delete [] sort_index;
}

size_t fastmarching::get_npesos()
{
	return n_pesos;
}

void fastmarching::print(ostream & os)
{
	for(size_t i = 0; i < n_vertices ; i++)
		os<<distances[i]<<endl;
}

void fastmarching::print_radio(ostream & os)
{
	for(size_t i = 0; i < n_vertices ; i++)
		if(distances[i] <= radio) os<<1<<" ";
		else os<<0<<" ";

	os<<endl;
}

distance_t fastmarching::operator[](size_t i)
{
	if(i < n_vertices)
		return distances[i];
	return -1;
}

size_t fastmarching::operator()(size_t i)
{
	if(i < n_pesos)
		return sort_index[i];
	return NIL;
}


void fastmarching::print_path(ostream & os)
{
	for(size_t i = 0; i < n_vertices ; i++)
		os<<paths[i].v<<endl;
}

void fastmarching::print_path(off & model, size_t x, ostream & os, size_t source)
{
	os<<model(x)<<endl;
	print_path(model, paths[x], os, distances[x] - *(model(x) - paths[x].v), source);
}

void fastmarching::print_path(off & model, path p, ostream & os, distance_t d, size_t source)
{
	os<<p.v<<endl;

	if(p.i == source)
	{
		//cout<<d<<" -> "<<distances[p.i]<<" ** "<<*(model(source) - p.v)<<" ***"<<endl;
		return;
	}
	if(p.i != NIL)
	{
		//cout<<d<<" -> "<<distances[p.i]<<" ** "<<*(model(source) - p.v)<<" ***"<<endl;
		print_path(model, paths[p.i], os, d - *(p.v - paths[p.i].v), source);
		return;
	}

	size_t vi;

	for(auto fi: model.get_faces(p.f[0]))
		for(auto fj: model.get_faces(p.f[1]))
			if(fi == fj)
			{
				if(model[fi][0] != p.f[2] && model[fi][1] != p.f[2] && model[fi][2] != p.f[2])
				{
					if(model[fi][0] != p.f[0] && model[fi][0] != p.f[1]) vi = model[fi][0];
					if(model[fi][1] != p.f[0] && model[fi][1] != p.f[1]) vi = model[fi][1];
					if(model[fi][2] != p.f[0] && model[fi][2] != p.f[1]) vi = model[fi][2];
				}
			}

	path r0, r1;
	face f0, f1;
	f0[2] = p.f[1];
	f1[2] = p.f[0];
	f0[0] = p.f[0];
	f0[1] = f1[0] = vi;
	f1[1] = p.f[1];

	distance_t p0 = update(path(f0, p.v, p.f[1]), r0, model);
	distance_t p1 = update(path(f1, p.v, p.f[0]), r1, model);

	//if(abs(p0 - d) < abs(p1 - d))
	if(p0 < p1)
	{
		//cout<<d<<" -> "<<p0<<" ** "<<*(model(source) - p.v)<<endl;
		print_path(model, r0, os, d - *(p.v - r0.v), source);
	}
	else
	{
		//cout<<d<<" -> "<<p1<<" ** "<<*(model(source) - p.v)<<endl;
		print_path(model, r1, os,  d - *(p.v - r1.v), source);
	}
}

/*
 * Color of vertices: black green red are  0 1 2  respectively
 */
void fastmarching::calcular_FM(off & model, vector<index_t> & source)
{
	for(index_t i = 0; i < model.get_nvertices(); i++)
	{
		distances[i] = INFINITY;
		color[i] = green;
	}

	size_t green_count = model.get_nvertices();
	priority_queue<pair<distance_t, size_t>,
			vector<pair<distance_t, size_t> >,
			greater<pair<distance_t, size_t> > > cola;

	path pr;
	distance_t p;
	size_t black_i;

	index_t c = 0;
	for(index_t s: source)
	{
		distances[s] = 0;
		cluster[s] = ++c;
		color[s] = red;
		cola.push(make_pair(distances[s], s));
	}

	while(green_count-- && !cola.empty())
	{
		while(!cola.empty() && color[cola.top().second] == black) cola.pop();
		if(cola.empty()) break;

		black_i = cola.top().second;
		color[black_i] = black;
		cola.pop();

		if(distances[black_i] > radio) break;

		sort_index[n_pesos++] = black_i;
	
		for(auto it: model.get_rings(black_i))
		{
			if(color[it] == green)
				color[it] = red;

			if(color[it] == red)
			{
				for(auto fi: model.get_faces(it))
				{
					p = update(path(model[fi], model(it), it), pr, model);
					
					if(p < distances[it])
					{
						distances[it] = p;
						cluster[it] = cluster[black_i];
						paths[it] = pr;
					}
				}
				if(distances[it] < INFINITY)
					cola.push(make_pair(distances[it], it));
			}
		}
	}
}

distance_t fastmarching::update(path p, path & r, off & model)
{
	mat X(3,2);
	vertex x0, x1, x2;	//x2 -> actualizar
	size_t ix0, ix1, ix2;

	ix2 = p.i;
	x2 = p.v;

	if(ix2 == p.f[2])
	{
		x0 = model(p.f[0]);
		x1 = model(p.f[1]);

		ix0 = p.f[0];
		ix1 = p.f[1];
	}

	if(ix2 == p.f[1])
	{
		x0 = model(p.f[0]);
		x1 = model(p.f[2]);

		ix0 = p.f[0];
		ix1 = p.f[2];
	}

	if(ix2 == p.f[0])
	{
		x0 = model(p.f[1]);
		x1 = model(p.f[2]);

		ix0 = p.f[1];
		ix1 = p.f[2];
	}

	vertex v[2];

	v[0] = x0 - x2;
	v[1] = x1 - x2;

	X(0, 0) = v[0][0];
	X(1, 0) = v[0][1];
	X(2, 0) = v[0][2];

	X(0, 1) = v[1][0];
	X(1, 1) = v[1][1];
	X(2, 1) = v[1][2];
	

	r.f[0] = ix0;
	r.f[1] = ix1;
	r.f[2] = ix2;
	r.v = p.v;
	r.i = NIL;

	if(spherical) return spherical_update(r, X, ix0, ix1);
	return planar_update(r, X, ix0, ix1);
}

distance_t fastmarching::planar_update(path & r, mat & X, size_t ix0, size_t ix1, bool debug)
{
	mat ones(2,1);
	ones.ones(2,1);

	mat Q;
	if(!inv_sympd(Q, X.t() * X))
		return INFINITY;
	
	mat t(2,1);

	t(0) = distances[ix0];
	t(1) = distances[ix1];

	distance_t p;
	mat delta = ones.t() * Q * t;
	distance_t dis = as_scalar(delta * delta - (ones.t() * Q * ones) * (as_scalar(t.t() * Q * t) - 1));

	if(dis >= 0)
	{
		p = delta(0) + sqrt(dis);
		p /= as_scalar(ones.t() * Q * ones);
		if(debug) cout<<"P: "<<p<<endl<<endl;
	}
	else p = INFINITY;

	mat n = X * Q * (t - p * ones);
	mat cond = Q * X.t() * n;

	mat x;
	if(t(0) == INFINITY || t(1) == INFINITY || dis < 0 || (cond(0) >= 0 || cond(1) >= 0))
	{
		if(debug) cout<<"-- dijkstra --"<<endl<<endl;
		distance_t p0 = distances[ix0] + norm(X.col(0));
		distance_t p1 = distances[ix1] + norm(X.col(1));

		if(p0 < p1)
		{
			r.i = ix0;
			x = X.col(0);
			p = p0;
		}
		else
		{
			r.i = ix1;
			x = X.col(1);
			p = p1;
		}
	}
	else
	{
		mat A(3,2);
		A.col(0) = -n;
		A.col(1) = X.col(1) - X.col(0);
		mat b = -X.col(0);
		mat l = solve(A, b);
		x = l(1) * A.col(1) + X.col(0);
	}

	r.v += vertex(x(0), x(1), x(2));

	return p;
}

distance_t fastmarching::spherical_update(path & r, mat & X, size_t ix0, size_t ix1)
{
	mat ones(2,1);
	ones.ones();

	mat Q;
	if(!inv_sympd(Q, X.t() * X))
		return INFINITY;
	
	mat s(2,1);

	s(0) = distances[ix0]*distances[ix0];
	s(1) = distances[ix1]*distances[ix1];

	mat q = s;

	q(0) -= as_scalar(X.col(0).t() * X.col(0));
	q(1) -= as_scalar(X.col(1).t() * X.col(1));

	distance_t p, dis;
	mat delta = ones.t()*Q*q;

	delta(0) += 2;
	dis = as_scalar(delta * delta - (ones.t() * Q * ones) * (q.t() * Q * q));

	if(dis >= 0)
	{
		p = delta(0) + sqrt(dis);
		p /= as_scalar(ones.t() * Q * ones);
	}
	else p = INFINITY;

	mat cond = Q * (p * ones - q);

	mat x;
	if(s(0) == INFINITY || s(1) == INFINITY || dis < 0 || cond(0) * cond(1) < 0)
	//if(s(0) == INFINITY || s(1) == INFINITY || dis < 0 || cond(0) >= 0 || cond(1) >= 0)
	{
		distance_t p0 = distances[ix0] + norm(X.col(0));
		distance_t p1 = distances[ix1] + norm(X.col(1));

		if(p0 < p1)
		{
			r.i = ix0;
			x = X.col(0);
			p = p0;
		}
		else
		{
			r.i = ix1;
			x = X.col(1);
			p = p1;
		}
	}
	else
	{
		mat n = solve(X.t(), p * ones - q);
		mat A(3,2);
		A.col(0) = -n;
		A.col(1) = X.col(1) - X.col(0);

		mat b = -X.col(0);

		mat l = solve(A, b);
		x = l(1)*A.col(1)+X.col(0);

		p = sqrt(p);
	}

	r.v += vertex(x(0), x(1), x(2));

	return p;
}

void fastmarching::normalize()
{
	distance_t max = distances[sort_index[n_pesos - 1]];
	for(size_t i = 0; i < n_pesos; i++)
		distances[sort_index[i]] /= max;
}
