#include <gproshan/laplacian/fairing_spectral.h>

#include <gproshan/laplacian/laplacian.h>


// geometry processing and shape analysis framework
namespace gproshan {


fairing_spectral::fairing_spectral(const size_t n_eigs_): n_eigs(n_eigs_) {}

void fairing_spectral::compute(che * mesh)
{
	delete [] vertices;
	vertices = new vertex[mesh->n_vertices];
	memcpy(vertices, &mesh->point(0), mesh->n_vertices * sizeof(vertex));

	arma::fmat X((float *) vertices, 3, mesh->n_vertices, false, true);

	arma::sp_fmat L, A;
	arma::fvec eigval;
	arma::fmat eigvec;

	n_eigs = eigs_laplacian(mesh, eigval, eigvec, L, A, n_eigs);
	if(!n_eigs) return;

	X = X * eigvec * eigvec.t();
}


} // namespace gproshan

