#include "laplacian/fairing_spectral.h"

#include "laplacian/laplacian.h"


// geometry processing and shape analysis framework
namespace gproshan {


fairing_spectral::fairing_spectral(const size_t & k_): fairing(), k(k_) {}

void fairing_spectral::compute(che * mesh)
{
	delete [] vertices;
	vertices = new vertex[mesh->n_vertices];
	memcpy(vertices, &mesh->gt(0), mesh->n_vertices * sizeof(vertex));

	a_mat X((real_t *) vertices, 3, mesh->n_vertices, false, true);

	a_sp_mat L, A;
	a_vec eigval;
	a_mat eigvec;

	k = eigs_laplacian(mesh, eigval, eigvec, L, A, k);

	X = X * eigvec * eigvec.t();
}


} // namespace gproshan

