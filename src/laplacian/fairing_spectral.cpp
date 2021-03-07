#include "laplacian/fairing_spectral.h"

#include "laplacian/laplacian.h"


// geometry processing and shape analysis framework
namespace gproshan {


fairing_spectral::fairing_spectral(const size_t & k_): fairing(), k(k_)
{
}

fairing_spectral::~fairing_spectral()
{

}

void fairing_spectral::compute(che * mesh)
{
	double time;

	positions = new vertex[mesh->n_vertices];

	a_mat X((real_t *) positions, 3, mesh->n_vertices, false, true);

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
		positions[v] = mesh->gt(v);

	a_sp_mat L, A;
	a_vec eigval;
	a_mat eigvec;

	TIC(time) k = eigs_laplacian(mesh, eigval, eigvec, L, A, k); TOC(time)
	gproshan_debug_var(time);

	X = X * eigvec * eigvec.t();
}


} // namespace gproshan

