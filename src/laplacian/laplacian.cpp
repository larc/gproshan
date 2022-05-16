#include "laplacian/laplacian.h"

using namespace std;
using namespace Eigen;


// geometry processing and shape analysis framework
namespace gproshan {


void laplacian(const che * mesh, a_sp_mat & L, a_sp_mat & A)
{
	size_t n_edges = mesh->n_edges;
	size_t n_vertices = mesh->n_vertices;

	arma::umat DI(2, 2 * n_edges);
	a_vec DV(2 * n_edges);

	arma::umat SI(2, n_edges);
	a_vec SV(n_edges);

	#pragma omp parallel for
	for(index_t e = 0; e < n_edges; ++e)
	{
		index_t i = e << 1;

		DI(0, i) = e;
		DI(1, i) = mesh->vt(mesh->et(e));
		DV(i) = -1;

		++i;

		DI(0, i) = e;
		DI(1, i) = mesh->vt(next(mesh->et(e)));
		DV(i) = 1;

		SI(0, e) = SI(1, e) = e;
		SV(e) = (mesh->cotan(mesh->et(e)) +
					mesh->cotan(mesh->ot_et(e))) / 2;
	}

	a_sp_mat D(DI, DV, n_edges, n_vertices);
	a_sp_mat S(SI, SV, n_edges, n_edges);

	L = D.t() * S * D;

	A.eye(n_vertices, n_vertices);

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		A(v, v) = mesh->area_vertex(v);
}

void laplacian(const che * mesh, sp_mat_e & L, sp_mat_e & A)
{
	gproshan_debug(LAPLACIAN);

	size_t n_edges = mesh->n_edges;
	size_t n_vertices = mesh->n_vertices;

	sp_mat_e D(n_edges, n_vertices);
	sp_mat_e S(n_edges, n_edges);

	D.reserve(VectorXi::Constant(n_edges,2));
	S.reserve(VectorXi::Constant(n_edges,1));

	for(index_t e = 0; e < n_edges; ++e)
	{
		D.insert(e, mesh->vt(mesh->et(e))) = 1;
		D.insert(e, mesh->vt(next(mesh->et(e)))) = -1;

		S.insert(e, e) = (mesh->cotan(mesh->et(e)) +
					mesh->cotan(mesh->ot_et(e))) / 2;
	}

	L = D.transpose() * S * D;

	A.reserve(VectorXi::Constant(n_vertices, 1));
	for(index_t v = 0; v < n_vertices; ++v)
		A.insert(v, v) = mesh->area_vertex(v);
		//A.insert(v, v) = 1.0 / sqrt(mesh->area_vertex(v));
}

size_t eigs_laplacian(const che * mesh, a_vec & eigval, a_mat & eigvec, a_sp_mat & L, a_sp_mat & A, const size_t & k)
{
	laplacian(mesh, L, A);

	string feigval = tmp_file_path(mesh->name_size() + ".eigval");
	string feigvec = tmp_file_path(mesh->name_size() + ".eigvec");

	if(!eigval.load(feigval) || !eigvec.load(feigvec) || eigval.n_elem < k)
	{
		if(!eigs_sym(eigval, eigvec, A * L * A, k, "sm"))
			return 0;

		eigval.save(feigval);
		eigvec.save(feigvec);
	}

	if(k < eigval.n_elem)
	{
		eigval = eigval.head(k);
		eigvec = eigvec.head_cols(k);
	}

	return eigval.n_elem;
}


} // namespace gproshan

