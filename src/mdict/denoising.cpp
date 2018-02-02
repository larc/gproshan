#include "denoising.h"

#include "sampling.h"
#include "d_dict_learning.h"

#include "che_poisson.h"
#include "che_fill_hole.h"
#include <cassert>
#include <set>

// mesh dictionary learning and sparse coding namespace
namespace mdict {

denoising::denoising(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _plot): dictionary(_mesh, _phi_basis, _m, _M, _f, _plot)
{
}

void denoising::execute()
{
	TIC(d_time) init_sampling(); TOC(d_time)
	debug(d_time)
	
	TIC(d_time) init_patches(); TOC(d_time)
	debug(d_time)
	
	TIC(d_time) learning(); TOC(d_time)
	debug(d_time)

	TIC(d_time) sparse_coding(); TOC(d_time)
	debug(d_time)

	d_message(mesh reconstruction...)
	assert(n_vertices == mesh->n_vertices());

	TIC(d_time)
	mesh_reconstruction(mesh, M, patches, patches_map, A, alpha);
	TOC(d_time)
	debug(d_time)
}

} // mdict

