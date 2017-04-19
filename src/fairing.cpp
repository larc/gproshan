#include "fairing.h"

fairing::fairing()
{
	positions = NULL;
}

fairing::~fairing()
{
	if(positions) delete [] positions;
}

void fairing::run(che * shape)
{
	compute(shape);
}

vertex * fairing::get_postions()
{
	return positions;
}
