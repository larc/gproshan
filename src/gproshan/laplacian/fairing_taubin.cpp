#include <gproshan/laplacian/fairing_taubin.h>

#include <gproshan/laplacian/laplacian.h>


// geometry processing and shape analysis framework
namespace gproshan {


fairing_taubin::fairing_taubin(const float step_): step(step_) {}

void fairing_taubin::compute(che * mesh)
{
	delete [] vertices;
	vertices = new vertex[mesh->n_vertices];
	memcpy(vertices, &mesh->point(0), mesh->n_vertices * sizeof(vertex));

	arma::fmat X((float *) vertices, 3, mesh->n_vertices, false, true);

	arma::sp_fmat L, A;
	laplacian(mesh, L, A);

	arma::fmat R = spsolve(A + step * L, A * X.t());

	X = R.t();
}


} // namespace gproshan

