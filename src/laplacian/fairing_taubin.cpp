#include "laplacian/fairing_taubin.h"

#include "laplacian/laplacian.h"


// geometry processing and shape analysis framework
namespace gproshan {


fairing_taubin::fairing_taubin(const real_t & step_): fairing(), step(step_)
{
}

fairing_taubin::~fairing_taubin()
{

}

void fairing_taubin::compute(che * shape)
{
	double time;
/*
	a_sp_mat_e Le, Ae;
	TIC(time)
	laplacian(shape, Le, Ae);
	TOC(time)
	cout<<"time laplacian: "<<time<<endl;
*/
	a_sp_mat L, A;

	gproshan_debug(compute laplacian);

	TIC(time) laplacian(shape, L, A); TOC(time)
	gproshan_debug_var(time);

	positions = new vertex[shape->n_vertices()];

	a_mat X((real_t *) positions, 3, shape->n_vertices(), false, true);

	#pragma omp parallel for
	for(index_t v = 0; v < shape->n_vertices(); v++)
		positions[v] = shape->gt(v);

	a_mat R;
	a_mat AX = A * X.t();
	a_sp_mat M = A + step * L;

	gproshan_debug(solve system);

	TIC(time) spsolve(R, M, AX); TOC(time)
	gproshan_debug_var(time);

	X = R.t();
}


} // namespace gproshan

