#include <gproshan/mesh/che_poisson.h>

#include <gproshan/mesh/toplesets.h>
#include <gproshan/laplacian/laplacian.h>

#include <armadillo>


// geometry processing and shape analysis framework
namespace gproshan {


void poisson(che * mesh, const size_t old_n_vertices, index_t k)
{
	if(!k) return;

	arma::fmat B(mesh->n_vertices, 3);
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		if(v < old_n_vertices)
		{
			B(v, 0) = mesh->point(v).x();
			B(v, 1) = mesh->point(v).y();
			B(v, 2) = mesh->point(v).z();
		}
		else B.row(v).zeros();
	}

	arma::sp_fmat L, A;
	laplacian(mesh, L, A);

	for(index_t i = 0; i < mesh->n_vertices; ++i)
		B.row(i) *= -1 / A(i,i);

	arma::sp_fmat M;
	float s = (k % 2) ? -1 : 1;

	if(k > 1) M = A * L;
	while(--k) L *= M;

	for(index_t v = 0; v < old_n_vertices; ++v)
	{
		B.col(0) -= mesh->point(v).x() * L.col(v);
		B.col(1) -= mesh->point(v).y() * L.col(v);
		B.col(2) -= mesh->point(v).z() * L.col(v);
	}

	L.shed_cols(0, old_n_vertices - 1);
	L.shed_rows(0, old_n_vertices - 1);
	A.shed_rows(0, old_n_vertices - 1);
	A.shed_cols(0, old_n_vertices - 1);
	B.shed_rows(0, old_n_vertices - 1);

	arma::fmat X;
	if(spsolve(X, s * L, s * B))
	for(index_t v = old_n_vertices; v < mesh->n_vertices; ++v)
	{
		mesh->point(v).x() = X(v - old_n_vertices, 0);
		mesh->point(v).y() = X(v - old_n_vertices, 1);
		mesh->point(v).z() = X(v - old_n_vertices, 2);
	}
}

void biharmonic_interp_2(arma::fmat & P, arma::fmat & H)
{
	size_t n = P.n_cols;
	float x;

	arma::fmat A(n, n);
	arma::fvec pi(2), pj(2);

	for(index_t i = 0; i < n; ++i)
	{
		pi(0) = P(0, i); pi(1) = P(1, i);
		for(index_t j = 0; j < n; ++j)
		{
			pj(0) = P(0, j); pj(1) = P(1, j);
			x = norm(pi - pj);
			A(i, j) = x * x * (log(x + 1e-18) - 1);
		}
	}

	arma::fmat alpha = solve(A, P.row(2).t());

	for(index_t i = 0; i < H.n_cols; ++i)
	{
		H(2, i) = 0;
		pi(0) = H(0, i); pi(1) = H(1, i);
		for(index_t j = 0; j < n; ++j)
		{
			pj(0) = P(0, j); pj(1) = P(1, j);
			x = norm(pi - pj);
			x *= x * (log(x + 1e-18) - 1);
			H(2, i) += alpha(j, 0) * x;
		}
	}
}

//fill one hole and fit with biharmonic_interp_2
void biharmonic_interp_2(che * mesh, const size_t old_n_vertices, const size_t n_vertices, const std::vector<index_t> & border_vertices, const index_t k)
{
	if(old_n_vertices == n_vertices) return;

	toplesets tps(mesh, border_vertices, k);

	const size_t n_border_vertices = tps.splits.back();

	std::vector<index_t> sub_mesh_hole;
	sub_mesh_hole.reserve(n_border_vertices);

	for(index_t b = 0; b < n_border_vertices; ++b)
		if(tps.sorted[b] < old_n_vertices)
			sub_mesh_hole.push_back(tps.sorted[b]);

	arma::fmat P(3, size(sub_mesh_hole));
	index_t i = 0;
	for(index_t & b: sub_mesh_hole)
	{
		P(0, i) = mesh->point(b).x();
		P(1, i) = mesh->point(b).y();
		P(2, i) = mesh->point(b).z();
		++i;
	}

	arma::fmat H(3, n_vertices - old_n_vertices);

	for(index_t i = 0, v = old_n_vertices; v < n_vertices; ++i, ++v)
	{
		H(0, i) = mesh->point(v).x();
		H(1, i) = mesh->point(v).y();
		H(2, i) = mesh->point(v).z();
	}

	arma::fvec avg = mean(H, 1);

	P.each_col() -= avg;
	H.each_col() -= avg;

	arma::fmat E;
	arma::fvec eval;
	eig_sym(eval, E, H * H.t());
	E.swap_cols(0,2);

	P = E.t() * P;
	H = E.t() * H;

	biharmonic_interp_2(P, H);

	H = E * H;
	H.each_col() += avg;

	mesh->update_vertices((vertex *) H.memptr(), H.n_cols, old_n_vertices);
}


} // namespace gproshan

