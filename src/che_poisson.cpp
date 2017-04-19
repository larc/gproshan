#include "che_poisson.h"
#include "laplacian.h"

#include <armadillo>

using namespace arma;

void poisson(che * mesh, const size_t & old_n_vertices, index_t k)
{
	if(!k) return;

	mat B(mesh->n_vertices(), 3);
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		if(v < old_n_vertices)
		{
			B(v, 0) = mesh->gt(v).x;
			B(v, 1) = mesh->gt(v).y;
			B(v, 2) = mesh->gt(v).z;
		}
		else B.row(v).zeros();
	}
	
	sp_mat L, A;
	laplacian(mesh, L, A);

	for(index_t i = 0; i < mesh->n_vertices(); i++)
		B.row(i) *= -1 / A(i,i);

	sp_mat M;
	vertex_t s = (k % 2) ? -1 : 1;
	
	if(k > 1) M = A * L; 
	while(--k) L *= M;

	for(index_t v = 0; v < old_n_vertices; v++)
	{
		B.col(0) -= mesh->gt(v).x * L.col(v);
		B.col(1) -= mesh->gt(v).y * L.col(v);
		B.col(2) -= mesh->gt(v).z * L.col(v);
	}	

	L.shed_cols(0, old_n_vertices - 1);
	L.shed_rows(0, old_n_vertices - 1);
	A.shed_rows(0, old_n_vertices - 1);
	A.shed_cols(0, old_n_vertices - 1);
	B.shed_rows(0, old_n_vertices - 1);

	mat X;
	if(spsolve(X, s * L, s * B))
	for(index_t v = old_n_vertices; v < mesh->n_vertices(); v++)
	{
		mesh->get_vertex(v).x = X(v - old_n_vertices, 0);
		mesh->get_vertex(v).y = X(v - old_n_vertices, 1);
		mesh->get_vertex(v).z = X(v - old_n_vertices, 2);
	}
}

void biharmonic_interp_2(mat & P, mat & H)
{
	size_t n = P.n_cols;
	vertex_t x;

	mat A(n, n);
	vec pi(2), pj(2);

	for(index_t i = 0; i < n; i++)
	{
		pi(0) = P(0, i); pi(1) = P(1, i);
		for(index_t j = 0; j < n; j++)
		{
			pj(0) = P(0, j); pj(1) = P(1, j);
			x = norm(pi - pj);
			A(i, j) = x * x * (log(x + 1e-18) - 1);
		}
	}
	
	mat alpha = solve(A, P.row(2).t());

	for(index_t i = 0; i < H.n_cols; i++)
	{
		H(2, i) = 0;
		pi(0) = H(0, i); pi(1) = H(1, i);
		for(index_t j = 0; j < n; j++)
		{
			pj(0) = P(0, j); pj(1) = P(1, j);
			x = norm(pi - pj);
			x *= x * (log(x + 1e-18) - 1);
			H(2, i) += alpha(j, 0) * x;
		}	
	}
}

//fill one hole and fit with biharmonic_interp_2
void biharmonic_interp_2(che * mesh, const size_t & old_n_vertices, const size_t & n_vertices, const vector<index_t> & border_vertices, const index_t & k)
{
	if(old_n_vertices == n_vertices) return;

	index_t * rings = new index_t[mesh->n_vertices()];
	index_t * sorted = new index_t[mesh->n_vertices()];
	vector<index_t> limites;
	mesh->sort_by_rings(rings, sorted, limites, border_vertices, k);

	const size_t n_border_vertices = limites.back();
	
	vector<index_t> sub_mesh_hole;
	sub_mesh_hole.reserve(n_border_vertices);

	for(index_t b = 0; b < n_border_vertices; b++)
		if(sorted[b] < old_n_vertices)
			sub_mesh_hole.push_back(sorted[b]);
	
	delete [] rings;
	delete [] sorted;
	
	mat P(3, sub_mesh_hole.size());
	index_t i = 0;
	for(index_t & b: sub_mesh_hole)
	{
		P(0, i) = mesh->gt(b).x;
		P(1, i) = mesh->gt(b).y;
		P(2, i) = mesh->gt(b).z;
		i++;
	}

	mat H(3, n_vertices - old_n_vertices);
	
	for(index_t i = 0, v = old_n_vertices; v < n_vertices; v++)
	{
		H(0, i) = mesh->gt(v).x;
		H(1, i) = mesh->gt(v).y;
		H(2, i) = mesh->gt(v).z;
		i++;
	}

	vec avg = mean(H, 1);

	P.each_col() -= avg;
	H.each_col() -= avg;

	mat E;
	vec eval;
	eig_sym(eval, E, H * H.t());
	E.swap_cols(0,2);
	
	P = E.t() * P;
	H = E.t() * H;

	biharmonic_interp_2(P, H);
	
	H = E * H;
	H.each_col() += avg;
	
	mesh->set_vertices((vertex *) H.memptr(), H.n_cols, old_n_vertices);
	
}

