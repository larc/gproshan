#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include <gproshan/mesh/che.h>

#include <armadillo>


// geometry processing and shape analysis framework
namespace gproshan {


template<class T>
void laplacian(const che * mesh, arma::SpMat<T> & L, arma::SpMat<T> & A)
{
	size_t n_edges = mesh->n_edges;
	size_t n_vertices = mesh->n_vertices;

	arma::umat DI(2, 2 * n_edges);
	arma::Col<T> DV(2 * n_edges);

	arma::umat SI(2, n_edges);
	arma::Col<T> SV(n_edges);

	#pragma omp parallel for
	for(index_t e = 0; e < n_edges; ++e)
	{
		index_t i = e << 1;

		DI(0, i) = e;
		DI(1, i) = mesh->edge_u(e);
		DV(i) = -1;

		++i;

		DI(0, i) = e;
		DI(1, i) = mesh->edge_v(e);
		DV(i) = 1;

		SI(0, e) = SI(1, e) = e;
		SV(e) = (mesh->cotan(mesh->edge_he_0(e)) + mesh->cotan(mesh->edge_he_1(e))) / 2;
	}

	arma::SpMat D(DI, DV, n_edges, n_vertices);
	arma::SpMat S(SI, SV, n_edges, n_edges);

	L = D.t() * S * D;

	A.eye(n_vertices, n_vertices);

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		A(v, v) = mesh->area_vertex(v);
}

template<class T>
size_t eigs_laplacian(const che * mesh, arma::Col<T> & eigval, arma::Mat<T> & eigvec,
										arma::SpMat<T> & L, arma::SpMat<T> & A, const size_t k)
{
	laplacian(mesh, L, A);

	std::string feigval = tmp_file_path(mesh->name_size() + ".eigval");
	std::string feigvec = tmp_file_path(mesh->name_size() + ".eigvec");

	if(!eigval.load(feigval) || !eigvec.load(feigvec) || eigval.n_elem < k)
	{
		if(!eigs_sym(eigval, eigvec, L, k, "sm"))
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

#endif // LAPLACIAN_H

