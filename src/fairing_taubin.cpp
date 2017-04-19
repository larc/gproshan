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

	mat X(shape->n_vertices(), 3);
	for(index_t v = 0; v < shape->n_vertices(); v++)
	{
		X(v, 0) = shape->gt(v).x;
		X(v, 1) = shape->gt(v).y;
		X(v, 2) = shape->gt(v).z;
	}

	mat AX = A * X;
	sp_mat M = A + step * L;
	
	d_message(Solve system...)
	TIC(time) spsolve(X, M, AX); TOC(time)
	debug(time)

	for(index_t v = 0; v < shape->n_vertices(); v++)
	{
		positions[v][0] = X(v, 0);
		positions[v][1] = X(v, 1);
		positions[v][2] = X(v, 2);
	}
}
