#include "synthesis.h"

// mesh dictionary learning and sparse coding namespace
namespace mdict {

synthesis::synthesis(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _plot): dictionary(_mesh, _phi_basis, _m, _M, _f, _plot)
{
}

void synthesis::execute()
{
	TIC(d_time) init_sampling(); TOC(d_time)
	debug(d_time)

	TIC(d_time) init_patches(); TOC(d_time)
	debug(d_time)

	//TIC(d_time) learning(); TOC(d_time)
	string name;
	d_message(Dictionary name:)
	cin>>name;
	string f_dict = "tmp/" + name + ".dict";
	debug(f_dict)
	if(!A.load(f_dict)) d_message(This dictionary does not exist Bye) return;
	debug(d_time)

	TIC(d_time) sparse_coding(); TOC(d_time)
	debug(d_time)

	TIC(d_time) mesh_reconstruction(); TOC(d_time)
	debug(d_time)
}

} // mdict

