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
	string f_dict;
	d_message(Dictionary name:)
	cin>>f_dict;
	//string f_dict = "tmp/" + mesh->name_size() + '_' + to_string(phi_basis->dim) + '_' + to_string(m) + ".dict";
	debug(f_dict)
	if(!A.load(f_dict)) d_message(This dictionary does not exist Bye)
	debug(d_time)

	TIC(d_time) sparse_coding(); TOC(d_time)
	debug(d_time)

	TIC(d_time) mesh_reconstruction(); TOC(d_time)
	debug(d_time)
}

} // mdict

