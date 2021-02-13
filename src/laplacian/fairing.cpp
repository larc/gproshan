#include "laplacian/fairing.h"


// geometry processing and shape analysis framework
namespace gproshan {


fairing::fairing()
{
	positions = nullptr;
}

fairing::~fairing()
{
	if(positions) delete [] positions;
}

void fairing::run(che * mesh)
{
	compute(mesh);
}

vertex * fairing::get_postions()
{
	return positions;
}


} // namespace gproshan

