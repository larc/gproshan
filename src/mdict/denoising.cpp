#include "denoising.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


denoising::denoising(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _learn, const bool & _plot): dictionary(_mesh, _phi_basis, _m, _M, _f, _learn, _plot)
{
}

distance_t denoising::execute()
{
	TIC(d_time) init_sampling(); TOC(d_time)
	gproshan_debug_var(d_time);

	TIC(d_time) init_patches(); TOC(d_time)
	gproshan_debug_var(d_time);

	TIC(d_time) learning(); TOC(d_time)
	gproshan_debug_var(d_time);

	TIC(d_time) sparse_coding(); TOC(d_time)
	gproshan_debug_var(d_time);
	draw_patches(128);
	draw_patches(12);
	phi_basis->plot_atoms(A);
	phi_basis->plot_basis();
	TIC(d_time) 
	distance_t error = mesh_reconstruction(); TOC(d_time)
	gproshan_debug_var(d_time);
	return error;
}


} // namespace gproshan::mdict

