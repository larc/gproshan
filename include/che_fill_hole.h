#ifndef CHE_FILL_HOLE_H
#define CHE_FILL_HOLE_H

#include "include.h"
#include "che.h"

#include <armadillo>

using namespace arma;

struct border_t
{
	angle_t theta;
	index_t v;

	border_t() = default;

	border_t(const vector<vec> & V, const index_t & _v, const array<index_t, 2> & neighbors, const bool & o):
	v(_v)
	{
		index_t p_v = neighbors[!o];
		index_t n_v = neighbors[o];
	
		if(p_v == NIL || n_v == NIL) 
		{
			theta = INFINITY;
			return;
		}

		vec a = V[p_v] - V[v];
		vec b = V[n_v] - V[v];
		a[2] = b[2] = 0;
		
		theta = atan2(b[1], b[0]) - atan2(a[1], a[0]);
		if(theta < 0) theta += 2 * M_PI;
	}

	vec new_vertex(const vector<vec> & V, angle_t div, const distance_t & lenght, const array<index_t, 2> & neighbors, const bool & o)
	{
		index_t p_v = neighbors[!o];
		
		vec a = V[p_v] - V[v];
		
		vec r(3);
		r[0] = cos(theta * div + atan2(a[1], a[0]));
		r[1] = sin(theta * div + atan2(a[1], a[0]));

		r = lenght * normalise(r);
		r += V[v];
		r[2] = 0;
		
		return r;
	}

	border_t(const vector<vec> & V, const index_t & _v, const array<index_t, 2> & neighbors, const bool & o, const vec & normal):
	v(_v)
	{
		index_t p_v = neighbors[!o];
		index_t n_v = neighbors[o];
		
		vec a = V[p_v] - V[v];
		vec b = V[n_v] - V[v];
	
		a -= dot(a, normal) * normal;
		b -= dot(a, normal) * normal;

		mat E(3,3);
		E.col(0) = normalise(a);
		E.col(1) = normalise(cross(normal, a));
		E.col(2) = normal;

		a = E.t() * a;
		b = E.t() * b;

		//theta = atan2( norm(cross(a,b)), dot(a,b));
		theta = atan2(b[1], b[0]) - atan2(a[1], a[0]);
		if(theta < 0) theta += 2 * M_PI;
	}

	vec new_vertex(const vector<vec> & V, angle_t div, const distance_t & lenght, const array<index_t, 2> & neighbors, const bool & o, const vec & normal)
	{
		index_t p_v = neighbors[!o];
		index_t n_v = neighbors[o];
		
		vec a = V[p_v] - V[v];
		vec b = V[n_v] - V[v];
		
		a -= dot(a, normal) * normal;
		b -= dot(a, normal) * normal;
		
		mat E(3,3);
		E.col(0) = normalise(a);
		E.col(1) = normalise(cross(normal, a));
		E.col(2) = normal;
		
		a = E.t() * a;
		b = E.t() * b;

		vec r(3);
		r[0] = cos(theta * div);
		r[1] = sin(theta * div);
		r[2] = 0;

		r = lenght * normalise(r);
		r = E * r;
		r += V[v];
		
		return r;
	}
};

bool operator<(const border_t & a, const border_t & b);

void poisson(che * mesh, const size_t & old_n_vertices);

vector<index_t> * fill_all_holes(che * mesh, const size_t & max_iter = 1000);

tuple<vector<index_t> *, che **> fill_all_holes_meshes(che * mesh, const size_t & max_iter = 1000);

che * fill_hole_front_angles_test(che * mesh, vector<index_t> & front_vertices, size_t p_iter, bool & is_grow);

che * fill_hole_front_angles_without_projection(che * mesh, vector<index_t> & front_vertices);

che * fill_hole_front_angles(vector<vertex> & vertices, const vertex_t & lenght, const vertex & normal, const size_t & max_iter, bool is_grow = false);

che * fill_hole_center_triangle(che * mesh, vector<index_t> & select_vertices, index_t index);

#endif //CHE_FILL_HOLE_H

