#include "fairing_taubin.h"
#include "laplacian.h"

fairing_taubin::fairing_taubin(matrix_t step_): fairing()
{
	 step = step_;
}

fairing_taubin::~fairing_taubin()
{

}

void fairing_taubin::compute(che * shape)
{
	double time;
/*
	sp_mat_e Le, Ae;
	TIC(time)
	laplacian(shape, Le, Ae);
	TOC(time)
	cout<<"time laplacian: "<<time<<endl;
*/
	sp_mat L, A;

	d_message(Compute laplacian...)
	TIC(time) laplacian(shape, L, A); TOC(time)
	debug(time)

	positions = new vertex[shape->n_vertices()];

	mat X((vertex_t *) positions, 3, shape->n_vertices(), false, true);

	#pragma omp parallel for
	for(index_t v = 0; v < shape->n_vertices(); v++)
		positions[v] = shape->gt(v);

	mat R;
	mat AX = A * X.t();
	sp_mat M = A + step * L;

	d_message(Solve system...)
	TIC(time) spsolve(R, M, AX); TOC(time)
	debug(time)

	X = R.t();
}

