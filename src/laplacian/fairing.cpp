#include "laplacian/fairing.h"


// geometry processing and shape analysis framework
namespace gproshan {


fairing::~fairing()
{
	delete [] vertices;
}

void fairing::run(che * mesh)
{
	compute(mesh);
}

const vertex * fairing::new_vertices()
{
	return vertices;
}


} // namespace gproshan

