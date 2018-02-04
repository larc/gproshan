#include "inpainting.h"

// mesh dictionary learning and sparse coding namespace
namespace mdict {

inpainting::inpainting(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _plot): dictionary(_mesh, _phi_basis, _m, _M, _f, _plot)
{
}

void inpainting::execute()
{

	n_vertices = mesh->n_vertices();
	delete [] fill_all_holes(mesh);
	
	TIC(d_time) poisson(mesh, n_vertices, 2); TOC(d_time)
	debug(d_time)

	debug(n_vertices)

	TIC(d_time) init_sampling(); TOC(d_time)
	debug(d_time)

// Here modify the patches function inside 

	TIC(d_time) init_patches(n_vertices); TOC(d_time)
	debug(d_time)
	
	TIC(d_time) learning(); TOC(d_time)
	debug(d_time)
	
	// Updating mesh n vertices after learning
	n_vertices = mesh->n_vertices();

	TIC(d_time) init_patches(); TOC(d_time)
	debug(d_time)

	TIC(d_time) sparse_coding(); TOC(d_time)
	debug(d_time)

	TIC(d_time) mesh_reconstruction(); TOC(d_time)
	debug(d_time)
}

} // mdict

