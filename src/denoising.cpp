#include "denoising.h"

#include "sampling.h"
#include "d_dict_learning.h"

#include "che_poisson.h"
#include "che_fill_hole.h"
#include <cassert>
#include <set>



denoising::~denoising()
{
}


void denoising::execute()
{
	d_message(sparse coding...)
	TIC(d_time)
	OMP_all_patches_ksvt(alpha, A, patches, M, L);
	TOC(d_time)

	d_message(mesh reconstruction...)
	assert(n_vertices == mesh->n_vertices());

	TIC(d_time)
	mesh_reconstruction(mesh, M, patches, patches_map, A, alpha);
	TOC(d_time)
	debug(d_time)
}
