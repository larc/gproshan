#include "inpainting.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


inpainting::inpainting(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _plot): dictionary(_mesh, _phi_basis, _m, _M, _f, _plot)
{
}

void inpainting::execute()
{
	// fill holes
	size_t threshold = mesh->n_vertices();
	delete [] fill_all_holes(mesh);
	TIC(d_time) poisson(mesh, threshold, 2); TOC(d_time)
	debug(d_time)

	// remove possible non manifold vertices
	mesh->remove_non_manifold_vertices();

	// sampling including new vertices
	TIC(d_time) init_sampling(); TOC(d_time)
	debug(d_time)

	// initializing patches with threshold
	TIC(d_time) init_patches(1, [&threshold](const index_t & i) -> bool { return i < threshold; }); TOC(d_time)
	debug(d_time)

	// learning only from valid patches
	TIC(d_time) learning(); TOC(d_time)
	debug(d_time)

	// including vertices out of threshold
	TIC(d_time) init_patches(0); TOC(d_time)
	debug(d_time)

	// Update new alphas, propagating the info towards the center
	TIC(d_time) update_alphas(alpha, threshold); TOC(d_time)
	debug(d_time)

	// sparse coding and reconstruction with all patches
	TIC(d_time) sparse_coding(); TOC(d_time)
	debug(d_time)

	TIC(d_time) mesh_reconstruction(); TOC(d_time)
	debug(d_time)
}


} // namespace gproshan::mdict

