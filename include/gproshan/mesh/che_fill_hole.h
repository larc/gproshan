#ifndef CHE_FILL_HOLE_H
#define CHE_FILL_HOLE_H

#include <gproshan/mesh/che.h>

#include <array>
#include <armadillo>


// geometry processing and shape analysis framework
namespace gproshan {


struct border_t
{
	float theta;
	index_t v;

	border_t() = default;

	border_t(const std::vector<arma::fvec> & V, const index_t _v, const std::array<index_t, 2> & neighbors, const bool o):
	v(_v)
	{
		index_t p_v = neighbors[!o];
		index_t n_v = neighbors[o];

		if(p_v == NIL || n_v == NIL)
		{
			theta = INFINITY;
			return;
		}

		arma::fvec a = V[p_v] - V[v];
		arma::fvec b = V[n_v] - V[v];
		a[2] = b[2] = 0;

		theta = atan2(b[1], b[0]) - atan2(a[1], a[0]);
		if(theta < 0) theta += 2 * M_PI;
	}

	arma::fvec new_vertex(const std::vector<arma::fvec> & V, float div, const float length, const std::array<index_t, 2> & neighbors, const bool o)
	{
		index_t p_v = neighbors[!o];
		index_t n_v = neighbors[o];

		arma::fvec a = V[p_v] - V[v];
		arma::fvec b = V[n_v] - V[v];

		a(2) = b(2) = 0;

		arma::fvec r = div * a + (1 - div) * b;

		r = length * normalise(r) + V[v];
		r(2) = 0;

		return r;
	}

	border_t(const std::vector<arma::fvec> & V, const index_t _v, const std::array<index_t, 2> & neighbors, const bool o, const arma::fvec & normal):
	v(_v)
	{
		index_t p_v = neighbors[!o];
		index_t n_v = neighbors[o];

		arma::fvec a = V[p_v] - V[v];
		arma::fvec b = V[n_v] - V[v];

		a -= dot(a, normal) * normal;
		b -= dot(a, normal) * normal;

		arma::fmat E(3,3);
		E.col(0) = normalise(a);
		E.col(1) = normalise(cross(normal, a));
		E.col(2) = normal;

		a = E.t() * a;
		b = E.t() * b;

		//theta = atan2( norm(cross(a,b)), dot(a,b));
		theta = atan2(b[1], b[0]) - atan2(a[1], a[0]);
		if(theta < 0) theta += 2 * M_PI;
	}

	arma::fvec new_vertex(const std::vector<arma::fvec> & V, float div, const float length, const std::array<index_t, 2> & neighbors, const bool o, const arma::fvec & normal)
	{
		index_t p_v = neighbors[!o];
		index_t n_v = neighbors[o];

		arma::fvec a = V[p_v] - V[v];
		arma::fvec b = V[n_v] - V[v];

		a -= dot(a, normal) * normal;
		b -= dot(a, normal) * normal;

		arma::fmat E(3,3);
		E.col(0) = normalise(a);
		E.col(1) = normalise(cross(normal, a));
		E.col(2) = normal;

		a = E.t() * a;
		b = E.t() * b;

		arma::fvec r(3);
		r[0] = cos(theta * div);
		r[1] = sin(theta * div);
		r[2] = 0;

		r = length * normalise(r);
		r = E * r;
		r += V[v];

		return r;
	}
};

bool operator<(const border_t & a, const border_t & b);

void poisson(che * mesh, const size_t old_n_vertices);

std::vector<index_t> * fill_all_holes(che * mesh, const size_t max_iter = 1000);

std::tuple<std::vector<index_t> *, che **> fill_all_holes_meshes(che * mesh, const size_t max_iter = 1000);

che * fill_hole_front_angles_test(che * mesh, std::vector<index_t> & front_vertices, size_t p_iter, bool & is_grow);

che * fill_hole_front_angles_without_projection(che * mesh, std::vector<index_t> & front_vertices);

che * fill_hole_front_angles(std::vector<vertex> & vertices, const float length, const vertex & normal, const size_t max_iter, bool is_grow = false);

che * fill_hole_center_triangle(che * mesh, std::vector<index_t> & select_vertices, index_t index);


} // namespace gproshan

#endif //CHE_FILL_HOLE_H

