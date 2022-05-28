#include <gproshan/laplacian/fairing_taubin.h>

#include <gproshan/laplacian/laplacian.h>


// geometry processing and shape analysis framework
namespace gproshan {


fairing_taubin::fairing_taubin(const real_t & step_): step(step_) {}

void fairing_taubin::compute(che * mesh)
{
	delete [] vertices;
	vertices = new vertex[mesh->n_vertices];
	memcpy(vertices, &mesh->gt(0), mesh->n_vertices * sizeof(vertex));

	a_mat X((real_t *) vertices, 3, mesh->n_vertices, false, true);

	a_sp_mat L, A;
	laplacian(mesh, L, A);

	a_mat R = spsolve(A + step * L, A * X.t());

	X = R.t();
}


} // namespace gproshan

