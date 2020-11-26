#include "laplacian/laplacian.h"

using namespace std;
using namespace Eigen;


// geometry processing and shape analysis framework
namespace gproshan {


void laplacian(che * mesh, a_sp_mat & L, a_sp_mat & A)
{
	size_t n_edges = mesh->n_edges();
	size_t n_vertices = mesh->n_vertices();

	arma::umat DI(2, 2 * n_edges);
	a_vec DV(2 * n_edges);

	arma::umat SI(2, n_edges);
	a_vec SV(n_edges);

	#pragma omp parallel for
	for(index_t e = 0; e < n_edges; e++)
	{
		index_t i = e << 1;

		DI(0, i) = e;
		DI(1, i) = mesh->vt(mesh->et(e));
		DV(i) = -1;

		i++;

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
	for(index_t v = 0; v < n_vertices; v++)
		A(v, v) = mesh->area_vertex(v);
}

void laplacian(che * mesh, sp_mat_e & L, sp_mat_e & A)
{
	gproshan_debug(LAPLACIAN);

	size_t n_edges = mesh->n_edges();
	size_t n_vertices = mesh->n_vertices();

	sp_mat_e D(n_edges, n_vertices);
	sp_mat_e S(n_edges, n_edges);

	D.reserve(VectorXi::Constant(n_edges,2));
	S.reserve(VectorXi::Constant(n_edges,1));

	for(index_t e = 0; e < n_edges; e++)
	{
		D.insert(e, mesh->vt(mesh->et(e))) = 1;
		D.insert(e, mesh->vt(next(mesh->et(e)))) = -1;

		S.insert(e, e) = (mesh->cotan(mesh->et(e)) +
					mesh->cotan(mesh->ot_et(e))) / 2;
	}

	L = D.transpose() * S * D;

	A.reserve(VectorXi::Constant(n_vertices, 1));
	for(index_t v = 0; v < n_vertices; v++)
		A.insert(v, v) = mesh->area_vertex(v);
}

size_t eigs_laplacian(a_vec & eigval, a_mat & eigvec, che * mesh, const a_sp_mat & L, const a_sp_mat & A, const size_t & K)
{
	gproshan_debug(LAPLACIAN);

	string feigval = tmp_file_path(mesh->name_size() + '_' + to_string(K) + ".L_eigval");
	string feigvec = tmp_file_path(mesh->name_size() + '_' + to_string(K) + ".L_eigvec");

	gproshan_debug_var(feigval);
	gproshan_debug_var(feigvec);

	if(!eigval.load(feigval) || !eigvec.load(feigvec))
	{
	//	a_sp_mat D = sqrt(A);
	//	D.for_each([](a_sp_mat::elem_type & val) { val = 1. / val; });

		if(!eigs_sym(eigval, eigvec, L, K, "sa"))
			return 0;
		
		eigval.save(feigval);
		eigvec.save(feigvec);
	}

	return eigval.n_elem;
}


} // namespace gproshan

