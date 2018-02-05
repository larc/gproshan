#include "fairing_spectral.h"
#include "laplacian.h"

fairing_spectral::fairing_spectral(size_t k_): fairing()
{
	k = k_;
}

fairing_spectral::~fairing_spectral()
{

}

void fairing_spectral::compute(che * shape)
{
	double time;

	sp_mat L, A;
	
	TIC(time) laplacian(shape, L, A); TOC(time)
	debug(time)

	positions = new vertex[shape->n_vertices()];

	mat X((vertex_t *) positions, 3, shape->n_vertices(), false, true);
	
	#pragma omp parallel for
	for(index_t v = 0; v < shape->n_vertices(); v++)
		positions[v] = shape->gt(v);

	vec eigval;
	mat eigvec;
	
	TIC(time) eigs_laplacian(eigval, eigvec, shape, L, k); TOC(time)
	debug(time)

	X = X * eigvec * eigvec.t();
}

